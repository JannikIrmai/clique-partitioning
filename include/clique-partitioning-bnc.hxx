#include <iostream>
#include <vector>
#include <iomanip>
#include <list>
#include <stack>

#include "gurobi_c++.h"

#include <kernighan-lin.hxx>
#include <property-maps.hxx>
#include <callback.hxx>
#include <time.hxx>

#include <rounding.hxx>

namespace CP {


/**
 * Node of the branch tree of a branch and cut algorithm.
 * It consists of an LP, pointers to its parent and children 
 * in the branch tree as well as some meta information like
 * i, j, value that describes which edge {i, j} was fixed
 * to which value or the depth of the node within the tree
 */
struct BranchNode
{
    // Constructs a new node with an empty LP
    BranchNode(GRBEnv env) : BranchNode(GRBModel(env)) {}

    // Construct a new branch node by copying an existing LP
    BranchNode(GRBModel m) : model(m) 
    {
        // set num-non-basic-iterations
        auto constraint_count = model.get(GRB_IntAttr_NumConstrs);
        GRBConstr* constraints = model.getConstrs();
        for (int i = 0; i < constraint_count; ++i) 
        {
            GRBConstr constraint = constraints[i];
            if (constraint.get(GRB_CharAttr_Sense) == '=')
                continue;
            num_non_basic_iterations_.push_back(0);
            constraints_.push_back(constraint);
        }
    }

    // Destructor destructs all children of the branch node and node itself
    ~BranchNode()
    {
        delete_children();
    }

    // Set the bound of the branch node and propagate this information to its parent:
    // The bound of the parent is the maximum of the bound of its two children.
    void set_bound(double new_bound)
    {
        assert (new_bound <= bound);
        if (parent)
        {
            double other_child_bound = value == 0 ? parent->child1->bound : parent->child0->bound;
            double new_parent_bound = std::max(other_child_bound, new_bound);
            if (new_parent_bound != parent->bound)
                parent->set_bound(new_parent_bound);
        }
        bound = new_bound;
    }

    // delete the children of the branch node
    void delete_children()
    {
        if (child0)
        {
            delete child0;
            child0 = nullptr;
        }
        if (child1)
        {
            delete child1;
            child1 = nullptr;
        }
    }

    // The number of (recursive) children of this node (including the node itself)
    size_t size()
    {
        size_t s = 1;
        if (child0)
            s += child0->size();
        if (child1)
            s += child1->size();
        return s;
    }

    // Creates two children of this branch node.
    // Both children are initialized with a copy of the LP of this branch node
    void create_children()
    {
        if (child0 || child1)
            throw std::runtime_error("BranchNode already has children!");

        child0 = new BranchNode(model);
        child0->parent = this;
        child0->depth = depth + 1;
        child0->bound = bound;

        child1 = new BranchNode(model);
        child1->parent = this;
        child1->depth = depth + 1;
        child1->bound = bound;
    }

    void print()
    {
        std::cout << "BranchNode at depth " << depth << " with i = " << i << ", j = " << j << ", and value = " << value << "\n";
    }

    size_t depth = 0;
    double bound;

    GRBModel model;

    size_t i;
    size_t j;
    size_t value;

    size_t min_stage_ = 0;

    BranchNode* parent = nullptr;
    BranchNode* child0 = nullptr;
    BranchNode* child1 = nullptr;

    bool lp_solved = false;

    // data for keeping track of which constraints are non basic
    std::list<size_t> num_non_basic_iterations_;
    std::list<GRBConstr> constraints_;

    // operator to compare two branch nodes. A branch node comes before another
    // if its bound is lower.
    bool operator<(const BranchNode& other) const
    {
        return bound < other.bound;
    }

    // prints the leaves of the branch tree routed at this branch node
    void print_leaves()
    {
        if (child0)
            child0->print_leaves();
        if (child1)
            child1->print_leaves();
        if (!child0 && !child1)
            std::cout << bound << " (" << depth << ") ";
    }
};

// method for comparing two branch nodes.
// Node n0 comes before n1 if its bound is smaller
auto branch_node_compare = [] (BranchNode* n0, BranchNode* n1) -> bool
{
    return n0->bound < n1->bound;
};

// Priority queue of branch nodes. 
// Top of the queue is the node with larges bound.
typedef std::priority_queue<BranchNode*, std::vector<BranchNode*>, decltype(branch_node_compare)> BranchNodeQueue;



/**
 * Class for managing the branch and cut algorithm for the clique partitioning problem.
 * This class maintains a branch and bound tree and (if branching is activated)
 * solves the clique partitioning problem by solving the linear relaxations of each branch node.
 * If the solution obtained by solving the LP of the current branch node is not integral, 
 * a fractional variable x_{i, j} is selected. Two children of the current branch node are created and
 * initialized with copies of the LP. The constraints x_{i, j} = 0 and x_{i, j} = 1 are added
 * to the respective LPs.
 * The constraints of the LP are modeled with lazily. That is, instead of enumerating
 * all inequalities upfront, only inequalities that are violated by the current solution to the LP
 * are added to the model. Moreover, inequalities that were not basic for several iterations 
 * are removed from the model in order to keep the model sparse.
 */
template<class EDGE_COST_MAP>
class BnC
{
public:
    typedef EdgePropertyMap<double> EPM;

    // public attributes
    // If a constraint is non-basic for more than max_iter_non_basic iterations, it is removed from the LP
    size_t max_iter_non_basic = std::numeric_limits<size_t>::max();
    // The higher the verbosity number is, the more information is printed to the console.
    size_t verbosity = 0;
    // The cutting plane algorithm (for solving an LP relaxation) is terminated if the relative decrease
    // of the objective value is less than tail_threshold for more than max_tail_length iterations
    double tail_threshold = 1.0;
    size_t max_tail_length = std::numeric_limits<size_t>::max();
    // if activate_branching is false, only the LP relaxation at the root node is solved.
    bool activate_branching = false;

    // Initialize an instance of the clique partitioning problem by specifying its costs
    BnC(const EDGE_COST_MAP& edge_costs) :
        BnC(edge_costs, GRBEnv())
    {}

    // Initialize an instance of the clique partitioning problem by specifying its costs
    BnC(const EDGE_COST_MAP& edge_costs, GRBEnv env) :
        n_(edge_costs.n()),
        separator_callback_(n_),
        edge_costs_(edge_costs),
        root_node_(env),
        branch_queue_(branch_node_compare),
        lp_solution_(n_),
        best_integer_solution_(n_)
    {
        root_node_.model.set(GRB_IntAttr_ModelSense, -1);  // maximize

        // add one variable for each edge {i, j} (and a copy for {j,i})
        vars_ = std::vector<std::vector<GRBVar>>(n_);
        for (size_t i = 0; i < n_; i++)
            vars_[i] = std::vector<GRBVar>(n_);

        for (size_t i = 0; i < n_; ++i)
        for (size_t j = i+1; j < n_; ++j)
        {
            vars_[i][j] = root_node_.model.addVar(0.0, 1.0, edge_costs(i, j), GRB_CONTINUOUS, 
                "x_"+std::to_string(i)+"_"+std::to_string(j));
            vars_[j][i] = vars_[i][j];
        }
        root_node_.model.set(GRB_IntParam_OutputFlag, 0);  // deactivate gurobi output
        root_node_.model.set(GRB_IntParam_Method, 1);  // use dual simplex
        current_node_ = &root_node_;
    }
    
    // This method adds some triangle inequalities to the model.
    // It can be beneficial to add those inequalities upfront instead of 
    // waiting until they are added lazily.
    void add_triangle_inequalities_based_on_negative_cost()
    {
        // for each edge ij with negative cost, add a triangle inequality of the
        // form x_ik + x_jk - x_ij <= 1. Pick k such that the minimum of the costs
        // ik and jk is maximal
        for (size_t i = 0; i < n_; ++i) 
        for (size_t j = i+1; j < n_; ++j)
        {
            if (edge_costs_(i, j) >= 0)
                continue;    
            size_t best_k;
            typename EDGE_COST_MAP::VALUE_TYPE max_min_cost = 0;
            for (size_t k = 0; k < n_; ++k)
            {
                if (k == i || k == j)
                    continue;
                auto min_cost = std::min(edge_costs_(i, k), edge_costs_(j, k));
                if (min_cost > max_min_cost)
                {
                    max_min_cost = min_cost;
                    best_k = k;
                }
            }
            if (max_min_cost > 0)
                add_constraint(vars_[i][best_k] + vars_[j][best_k] - vars_[i][j] <= 1);
        }
    }

    // This method adds some triangle inequalities to the model.
    // It can be beneficial to add those inequalities upfront instead of 
    // waiting until they are added lazily.
    void add_triangle_inequalities_based_on_positive_cost()
    {
        // for each edge ij with positive cost, add a triangle inequality of the
        // form x_ij + x_jk - x_jk <= 1. Pick k such that the minimum of the costs
        // ik and -jk is maximal
        for (size_t i = 0; i < n_; ++i) 
        for (size_t j = i+1; j < n_; ++j)
        {
            if (edge_costs_(i, j) <= 0)
                continue;    
            size_t best_k;
            bool best_k_flag;
            typename EDGE_COST_MAP::VALUE_TYPE max_min_cost = 0;
            for (size_t k = 0; k < n_; ++k)
            {
                if (k == i || k == j)
                    continue;
                auto min_cost = std::min(edge_costs_(i, k), -edge_costs_(j, k));
                if (min_cost > max_min_cost)
                {
                    max_min_cost = min_cost;
                    best_k = k;
                    best_k_flag = true;
                }
                min_cost = std::min(-edge_costs_(i, k), edge_costs_(j, k));
                if (min_cost > max_min_cost)
                {
                    max_min_cost = min_cost;
                    best_k = k;
                    best_k_flag = false;
                }
            }
            if (max_min_cost > 0)
            {
                if (best_k_flag)
                    add_constraint(vars_[i][j] + vars_[i][best_k] - vars_[j][best_k] <= 1);
                else
                    add_constraint(vars_[i][j] - vars_[i][best_k] + vars_[j][best_k] <= 1);
            }
        }
    }

    // This function optimizes the model.
    void optimize(double time_limit = std::numeric_limits<double>::max(), double best_objective = 0)
    {   
        optimization_start_time_ = Time::now();
        // compute trivial bound
        double trivial_bound = 0;
        for (size_t i = 0; i < n_; ++i) 
        for (size_t j = i+1; j < n_; ++j)
        {
            if (edge_costs_(i, j) > 0)
                trivial_bound += edge_costs_(i, j);
        }
        root_node_.bound = trivial_bound;

        // if branching is activated, an integer feasible solution is computed by means
        // of a local search heuristic (kernighan lin based greedy moving)
        best_objective_ = best_objective;
        if (activate_branching)
        {
            std::vector<size_t> node_labels(n_);
            std::iota(node_labels.begin(), node_labels.end(), 0);
            double kl_obj = kernighanLin(edge_costs_, node_labels);
            if (kl_obj > best_objective_)
            {
                best_objective_ = kl_obj;
                for (size_t i = 0; i < n_; ++i)
                for (size_t j = i+1; j < n_; ++j)
                    best_integer_solution_(i, j) = node_labels[i] == node_labels[j];
            }
        }

        // logging
        stage_ = 0;
        explored_node_count_ = 0;
        lp_time_ = 0;
        log_();

        time_limit_ = time_limit;

        print_status_bar_(std::cout);
        
        // reset the branch tree
        root_node_.delete_children();
        if (branch_queue_.size() != 0)
            throw std::runtime_error("branch_queue_ is not empty!");
        branch_queue_.push(&root_node_);

        while (true)
        {
            if (branch_queue_.size() == 0)
                throw std::runtime_error("Branch queue is empty.");
            // get the branch node with the largest bound
            current_node_ = branch_queue_.top();
            branch_queue_.pop();
            // assert that the bound is at least as strong as the bound of the root node
            assert (std::abs(current_node_->bound < root_node_.bound + 1e-3));
            // If the bound is less than one greater than the best objective, 
            // then the best objective is optimal
            if (current_node_->bound < best_objective_ + 1 - EPSILON)
                break;  // found optimal solution
            
            // if the lp has already been solved, branch
            if (current_node_->lp_solved)
            {
                if (!activate_branching)
                    break;
                branch_();
                continue;
            }
            // otherwise, solve the lp of the current node
            solve_current_node_lp_();
            if(current_node_->model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                break;

            // check if solution is integer
            if (num_integer() == n_ * (n_ - 1) / 2)
            {
                // update the best integer objective
                if (current_node_->bound > best_objective_)
                {
                    best_objective_ = current_node_->bound;
                    for (size_t i = 0; i < n_; ++i)
                    for (size_t j = i+1; j < n_; ++j)
                        best_integer_solution_(i, j) = lp_solution_(i, j);
                }
            }
            else if (current_node_->bound >= best_objective_ + 1 - EPSILON)
                branch_queue_.push(current_node_);  // If solution remains fractional, add it back to the queue for branching

            if (!activate_branching)
                break;

            // round the fractional solution in order to find a new best integer feasible solution
            if (current_node_->lp_solved)
            {
                std::vector<size_t> node_labels(n_);
                double obj = round_kl(lp_solution_, edge_costs_, 11, node_labels);
                if (obj > best_objective_)
                {
                    best_objective_ = obj;
                    for (size_t i = 0; i < n_; ++i)
                    for (size_t j = i+1; j < n_; ++j)
                        best_integer_solution_(i, j) = node_labels[i] == node_labels[j];
                }
            }

            if (root_node_.bound < best_objective_ + 1 - EPSILON)
            {
                break;  // found and proved optimal solution
                
            }
            if (((float_time_point)Time::now() - optimization_start_time_).count() > time_limit_)
                break;  // hit time limit
        }
        log_();
        // print the final status
        print_status_(std::cout);
        if (verbosity == 1)
            std::cout << "\n";
        // clear the branch tree
        root_node_.delete_children();
        current_node_ = &root_node_;
    }

    // return a reference to the separator callback (e.g. for adding additional separation routines)
    SeparatorCallback<EPM>& separator_callback()
    {
        return separator_callback_;
    }

    // return the global bound (i.e. bound at root node)
    double bound()
    {
        return root_node_.bound;
    }

    // return the best integer feasible objective
    double best_objective()
    {
        return best_objective_;
    }

    // return the number of integer variables in the current LP
    size_t num_integer()
    {
        size_t count = 0;
        for (size_t i = 0; i < n_; ++i)
        for (size_t j = i+1; j < n_; ++j)
        {
            double x = edge_value_(i, j);
            count += std::abs(x - int(x + 0.5)) < 1e-3;
        }
        return count;
    }

    // Print the variable values of the solution
    template<class OUT_STREAM>
    void print_solution(OUT_STREAM& out, bool integer=true)
    {
        out << n_ << "\n";
        for (size_t i = 0; i < n_; ++i)
        {
            for (size_t j = i+1; j < n_; ++j)
            {   
                double val = get_solution(i, j, integer);
                // round to 0 or 1 if close
                if (val < EPSILON)
                    val = 0;
                else if (val > 1 - EPSILON)
                    val = 1;
                out << val << " ";
            }
            out << "\n";
        }
    }

    double get_solution(size_t i, size_t j, bool integer)
    {
        return integer ? best_integer_solution_(i, j) : lp_solution_(i, j);
    }

    const std::vector<double>& get_log_time() const
    {
        return log_elapsed_time_;
    }

    const std::vector<double>& get_log_global_bound() const
    {
        return log_global_bound_;
    }

    const std::vector<double>& get_log_current_node_bound() const
    {
        return log_current_node_bound_;
    }

    const std::vector<double>& get_log_best_objective() const
    {
        return log_best_objective_;
    }

    const std::vector<double>& get_log_stage() const
    {
        return log_stages_;
    }

    const std::vector<double>& get_log_explored_node_count() const
    {
        return log_explored_node_count_;
    }

    const std::vector<double>& get_log_lp_time() const
    {
        return log_lp_time_;
    }

private:
    size_t n_;  // number of nodes of the instance
    EDGE_COST_MAP edge_costs_;  // edge costs of the instance

    std::vector<std::vector<GRBVar>> vars_;  // variables of the gurobi model
    
    EPM lp_solution_;  // variable values of the LP that has been solved last
    EPM best_integer_solution_;  // variable value of best integer feasible solution


    SeparatorCallback<EPM> separator_callback_;  // separator callback to separate fractional solution
    
    // time limits
    float_time_point optimization_start_time_;
    double time_limit_;
    double lp_time_;

    size_t num_separation_calls_ = 0;  // number of times the separator was called
    double best_objective_ = 0;  // maximum objective value of an integer feasible solution

    size_t tail_counter_;  // number of iterations in which the LP did not improve more than a certain threshold
    size_t stage_;  // separator stage of the last call back

    // vectors for logging different numerical value (time, bounds, etc.)
    std::vector<double> log_elapsed_time_;
    std::vector<double> log_current_node_bound_;
    std::vector<double> log_global_bound_;
    std::vector<double> log_best_objective_;
    std::vector<double> log_stages_;
    std::vector<double> log_lp_time_;
    std::vector<double> log_explored_node_count_;

    // representation of the branch tree
    BranchNode root_node_;
    BranchNode* current_node_;
    BranchNodeQueue branch_queue_;
    size_t explored_node_count_;

    // status bar of printout during branch and cut algorithm
    template<class OUT_STREAM>
    void print_status_bar_(OUT_STREAM& out)
    {
        if (verbosity == 0)
            return;  // do not print anything if verbosity is 0
        std::cout 
            << " Iter "
            << "EXPND "
            << "OPNND "
            << "DEPTH "
            << "    TIME "
            << " LP-TIME "
            << "   OBJBST "
            << "   OBJBND "
            << "   NODBND "
            << "  %I "
            << "#Constr ";
        for (size_t i = 0; i < separator_callback_.num_separators(); ++i)
            std::cout << std::setw(20) << separator_callback_.name(i) << " ";
        std::cout << "\n";
    }

    // print the numerical values below the status bar
    template<class OUT_STREAM>
    void print_status_(OUT_STREAM& out)
    {
        if (verbosity == 0)
            return;
        std::cout << "\r";  // overwrite the current status
        std::cout << std::setfill(' ') << std::setw(5) << num_separation_calls_ << " ";
        std::cout << std::setfill(' ') << std::setw(5) << explored_node_count_ << " ";
        std::cout << std::setfill(' ') << std::setw(5) << branch_queue_.size() << " ";
        std::cout << std::setfill(' ') << std::setw(5) << current_node_->depth << " ";
        double elapsed_time = ((float_time_point)Time::now() - optimization_start_time_).count();
        std::cout << std::setw(8) << std::setprecision(5) << elapsed_time << " ";
        std::cout << std::setw(8) << std::setprecision(5) << lp_time_ << " ";
        std::cout << std::setw(9) << std::setprecision(8) << best_objective() << " ";
        std::cout << std::setw(9) << std::setprecision(8) << bound() << " ";
        std::cout << std::setw(9) << std::setprecision(8) << current_node_->bound << " ";
        double frac_int = (double)num_integer() / (n_ * (n_ - 1) / 2);
        std::cout << std::setw(4) << std::setprecision(3) << 100 * frac_int << " ";
        std::cout << std::setw(7) << current_node_->model.get(GRB_IntAttr_NumConstrs) << " ";
        for (size_t i = 0; i < separator_callback_.num_separators(); ++i)
        {

            std::cout 
                << std::setw(4) << separator_callback_.num_calls(i) << " "
                << std::setw(7) << separator_callback_.num_inequalities(i) << " " 
                << std::setw(7) << std::setprecision(4) << separator_callback_.time(i) << " ";
        }
        if (verbosity >= 2)
            std::cout << "\n"; // if verbosity is large, print a new line instead of overwriting the previous line
        else 
            std::cout << "    " << std::flush;

    }

    // add a constraint to the LP of the current node of the branch tree
    void add_constraint(const GRBTempConstr& constr)
    {
        current_node_->constraints_.emplace_back(current_node_->model.addConstr(constr));
        current_node_->num_non_basic_iterations_.push_back(0);
    }

    // add a constraint to the LP of the current node of the branch tree
    void add_constraint(const Inequality<int>& inequality)
    {
        // convert the Inequality to a gurobi constraint
        GRBLinExpr expr;
        for (size_t i = 0; i < inequality.edges().size(); ++i)
        {
            size_t u = inequality.edges()[i][0];
            size_t v = inequality.edges()[i][1];
            expr += inequality.coefficients()[i] * vars_[u][v];
        }
        add_constraint(expr <= inequality.rhs());
    }

    // return the variable value x_{i, j} of the current LP
    double edge_value_(size_t i, size_t j)
    {
        GRBVar* var = &vars_[i][j];
        double* val;
        val = current_node_->model.get(GRB_DoubleAttr_X, var, 1);
        return val[0];
    }

    // separate the current lp solution and return the number of added constraints
    size_t separate_()
    {   
        ++num_separation_calls_;
        // call separator_callback_ to get violated inequalities
        auto inequalities = separator_callback_(lp_solution_, stage_, current_node_->min_stage_);
        // add those inequalities to the model
        for (const auto& inequality : inequalities)
            add_constraint(inequality);
        // return number of added inequalities
        return inequalities.size();
    }

    // remove constraints from the model that have been non-basic for a certain number of iterations
    size_t remove_non_basic_constraints()
    {
        assert (current_node_->constraints_.size() + current_node_->depth == current_node_->model.get(GRB_IntAttr_NumConstrs));
        assert (current_node_->num_non_basic_iterations_.size() == current_node_->constraints_.size()); 
        size_t num_removed = 0;  // count how many constraints are removed
        // iterate over constraints in current model
        auto it_constr = current_node_->constraints_.cbegin();
        auto it_count = current_node_->num_non_basic_iterations_.begin();
        while ( it_constr != current_node_->constraints_.cend() ) 
        {   
            // increase counter if constraint is non-basic
            if (it_constr->get(GRB_IntAttr_CBasis) == 0)
                ++(*it_count);
            // otherwise, reset the counter
            else
                *it_count = 0;
            
            // if the counter reached maximum value, remove constraint
            if (*it_count >= max_iter_non_basic)
            {
                current_node_->model.remove(*it_constr);
                num_removed++;
                it_constr = current_node_->constraints_.erase(it_constr);
                it_count = current_node_->num_non_basic_iterations_.erase(it_count);
            }
            else
            {
                ++it_constr;
                ++it_count;
            }
        }
        return num_removed;
    }

    // preform branching step of branch and cut algorithm:
    //  1. select an edge variable x_{i, j} for branching
    //  2. create two children of the current branch node and
    //      initialize them with copies of the current LP
    //  3. Add x_{i, j} = 0 and x_{i, j} = 1 to the two copies, respectively
    void branch_()
    {
        // assert that the current node has a fractional variable
        assert (num_integer() < n_ * (n_-1) / 2);
        // select a variable to branch on
        size_t i;
        size_t j; 
        select_branching_variable_(i, j);
        // if verbosity is very high, print variable
        if (verbosity >= 100)
            std::cout << "branching on " << i << " " << j << "\n";

        // add children to the current node
        current_node_->create_children();

        current_node_->child0->i = i;
        current_node_->child0->j = j;
        current_node_->child0->value = 0;
        current_node_->child0->model.addConstr(vars_[i][j] == 0);

        current_node_->child1->i = i;
        current_node_->child1->j = j;
        current_node_->child1->value = 1;
        current_node_->child1->model.addConstr(vars_[i][j] == 1);

        branch_queue_.push(current_node_->child0);
        branch_queue_.push(current_node_->child1);
    }

    // select for branching the edge {i,j} that maximizes
    //      x_{i,j} * (1-x_{i,j}) * c_{i, j}
    // However, if the maximal value is 0, then select the variable
    // that is closes to 0.5.
    void select_branching_variable_(size_t& i, size_t& j)
    {
        double largest_value = 0;
        size_t most_fractional_i;
        size_t most_fractional_j;
        double least_dist_to_one_half = 1;
        for (size_t a = 0; a < n_; ++a)
        for (size_t b = a+1; b < n_; ++b)
        {
            double x_ab = edge_value_(a, b);
            if (x_ab <= 0.001 || x_ab >= 0.999)
                continue; // x_ab is integral
            double value = x_ab * (1-x_ab) * edge_costs_(a, b);
            if (value > largest_value)
            {
                i = a;
                j = b;
                largest_value = value;
            }
            double dist_to_one_half = std::abs(0.5 - x_ab);
            if (dist_to_one_half < least_dist_to_one_half)
            {
                least_dist_to_one_half = dist_to_one_half;
                most_fractional_i = a;
                most_fractional_j = b;
            }
        }

        if (largest_value < EPSILON)
        {
            i = most_fractional_i;
            j = most_fractional_j;
        }
    }

    // Solve the LP of the current branch node
    // This is done by iteratively solving the LP, adding new violated constraints
    // and removing constraints that have been non-basic for a certain number of iterations
    void solve_current_node_lp_()
    {
        current_node_->model.update();
        tail_counter_ = 0;
        stage_ = 0;
        log_();
        while (true)
        {
            // solve the LP within time limit
            float_time_point lp_start = Time::now();
            double elapsed_time = (lp_start - optimization_start_time_).count();
            if (elapsed_time > time_limit_)
                break;
            current_node_->model.set(GRB_DoubleParam_TimeLimit, time_limit_ - elapsed_time);
            current_node_->model.optimize();
            float_time_point lp_end = Time::now();
            lp_time_ += (lp_end - lp_start).count();
            // check LP status
            if (current_node_->model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
                throw std::runtime_error("Model is infeasible!");
            if (current_node_->model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                return;
            // extract the solution from the LP
            for (size_t i = 0; i < n_; ++i)
            for (size_t j = i+1; j < n_; ++j)
                lp_solution_(i, j) = edge_value_(i, j);
            // extract the objective value
            double new_objective = current_node_->model.get(GRB_DoubleAttr_ObjVal);
            double old_objective = current_node_->bound;
            current_node_->set_bound(new_objective);

            // logging
            print_status_(std::cout);
            log_();

            // update the minimum stage of the separator for this node
            current_node_->min_stage_ = std::max(stage_, current_node_->min_stage_);

            // check if tailing of condition is met
            if (new_objective >= tail_threshold * old_objective)
            {
                ++tail_counter_;
                if (tail_counter_ >= max_tail_length)
                {
                    // If the minimum callback stage can be increase, increase the stage
                    if (current_node_->min_stage_ + 1 < separator_callback_.num_stages())
                    {
                        ++current_node_->min_stage_;
                        tail_counter_ = 0;
                    }
                    // otherwise terminate
                    else
                        break;
                }
            }
            else
            {
                tail_counter_ = 0;
            }

            // if objective increased, remove non-binding cuts
            if (new_objective < old_objective - EPSILON)
            {
                size_t num_removed = remove_non_basic_constraints();
            }

            // if the objective is less than one greater than the best integer feasible objective
            // then, this branch node need not to be considered any further
            if (activate_branching && current_node_->bound < best_objective_ + 1 - EPSILON)
                break;

            // call the separator callback
            size_t n = separate_();
            // if no additional inequalities are added to the model, this branch node is solved
            if (n == 0)
                break;
        }
        current_node_->lp_solved = true;
        ++explored_node_count_;
    }

    // log various numerical values
    void log_()
    {
        log_elapsed_time_.push_back(((float_time_point)Time::now() - optimization_start_time_).count());
        log_stages_.push_back(stage_);
        log_global_bound_.push_back(root_node_.bound);
        log_current_node_bound_.push_back(current_node_->bound);
        log_best_objective_.push_back(best_objective_);
        log_lp_time_.push_back(lp_time_);
        log_explored_node_count_.push_back(explored_node_count_);
    }
};

} // namespace CP

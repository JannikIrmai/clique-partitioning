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


// branching stuff
struct BranchNode
{
    BranchNode(GRBEnv env) : BranchNode(GRBModel(env)) {}
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
    ~BranchNode()
    {
        delete_children();
    }

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

    size_t size()
    {
        size_t s = 1;
        if (child0)
            s += child0->size();
        if (child1)
            s += child1->size();
        return s;
    }


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

    bool operator<(const BranchNode& other) const
    {
        return bound < other.bound;
    }

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

auto branch_node_compare = [] (BranchNode* n0, BranchNode* n1) -> bool
{
    return n0->bound < n1->bound;
};

typedef std::priority_queue<BranchNode*, std::vector<BranchNode*>, decltype(branch_node_compare)> BranchNodeQueue;





template<class EDGE_COST_MAP>
class BnC
{
public:
    typedef std::vector<double> VECTOR;
    typedef EdgePropertyMap<VECTOR> EPM;

    // public attributes
    size_t max_iter_non_basic = std::numeric_limits<size_t>::max();
    size_t verbosity = 0;
    double tail_threshold = 1.0;
    size_t max_tail_length = std::numeric_limits<size_t>::max();
    bool activate_branching = false;

    BnC(EDGE_COST_MAP edge_costs) :
        BnC(edge_costs, GRBEnv())
    {}

    BnC(EDGE_COST_MAP edge_costs, GRBEnv env) :
        n_(edge_costs.n()),
        callback_(n_),
        edge_costs_(edge_costs),
        root_node_(env),
        branch_queue_(branch_node_compare)
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
        root_node_.model.set(GRB_IntParam_OutputFlag, 0);
        root_node_.model.set(GRB_IntParam_Method, 1);  // use dual simplex
        current_node_ = &root_node_;
    }
    
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

        best_objective_ = best_objective;
        if (activate_branching)
        {
            std::vector<size_t> node_labels(n_);
            std::iota(node_labels.begin(), node_labels.end(), 0);
            double kl_obj = kernighanLin(edge_costs_, node_labels);
            if (kl_obj > best_objective_)
                best_objective_ = kl_obj;
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
            //std::cout << "queue.size() = " << branch_queue_.size() << "\n";
            // get the branch node with the largest bound
            current_node_ = branch_queue_.top();
            branch_queue_.pop();
            if (std::abs(current_node_->bound - root_node_.bound) > 1e-3)
            {
                throw std::runtime_error("This should not have happened!");
            }
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
                if (current_node_->bound > best_objective_)
                    best_objective_ = current_node_->bound;
            }
            else if (current_node_->bound >= best_objective_ + 1 - EPSILON)
                branch_queue_.push(current_node_);  // If solution remains fractional, add it back to the queue for branching

            if (!activate_branching)
                break;

            // round the fractional solution
            if (current_node_->lp_solved)
            {
                VECTOR edge_value_vector(n_ * (n_-1) / 2);
                EPM edge_value_map(n_, edge_value_vector);
                for (size_t i = 0; i < n_; ++i)
                for (size_t j = i+1; j < n_; ++j)
                    edge_value_map(i, j) = edge_value(i, j);
                double obj = round_kl(edge_value_map, edge_costs_, 11);
                if (obj > best_objective_)
                    best_objective_ = obj;
            }

            if (root_node_.bound < best_objective_ + 1 - EPSILON)
                break;  // found and proved optimal solution
            if (((float_time_point)Time::now() - optimization_start_time_).count() > time_limit_)
                break;  // hit time limit
        }
        log_();
        // print the final status
        print_status_(std::cout);
        if (verbosity == 1)
            std::cout << "\n";
        root_node_.delete_children();
    }

    Callback<EPM>& separator_callback()
    {
        return callback_;
    }

    double bound()
    {
        return root_node_.bound;
    }

    double best_objective()
    {
        return best_objective_;
    }

    double edge_value(size_t i, size_t j)
    {
        GRBVar* var = &vars_[i][j];
        double* val;
        val = current_node_->model.get(GRB_DoubleAttr_X, var, 1);
        return val[0];
    }

    size_t num_integer()
    {
        size_t count = 0;
        for (size_t i = 0; i < n_; ++i)
        for (size_t j = i+1; j < n_; ++j)
        {
            double x = edge_value(i, j);
            count += std::abs(x - int(x + 0.5)) < 1e-3;
        }
        return count;
    }

    template<class OUT_STREAM>
    void print_solution(OUT_STREAM& out)
    {
        out << n_ << "\n";
        for (size_t i = 0; i < n_; ++i)
        {
            for (size_t j = i+1; j < n_; ++j)
            {
                out << edge_value(i, j) << " ";
            }
            out << "\n";
        }
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
    size_t n_;
    EDGE_COST_MAP edge_costs_;

    std::vector<std::vector<GRBVar>> vars_;

    Callback<EPM> callback_;
    
    float_time_point optimization_start_time_;
    double time_limit_;
    double lp_time_;
    size_t num_separation_calls_ = 0;
    double best_objective_ = 0;

    size_t tail_counter_;
    size_t stage_;  // separator stage of the last call back
    size_t min_stage_;

    std::ofstream log_file_;

    // logging stuff
    std::vector<double> log_elapsed_time_;
    std::vector<double> log_current_node_bound_;
    std::vector<double> log_global_bound_;
    std::vector<double> log_best_objective_;
    std::vector<double> log_stages_;
    std::vector<double> log_lp_time_;
    std::vector<double> log_explored_node_count_;

    // branching stuff
    BranchNode root_node_;
    BranchNode* current_node_;
    BranchNodeQueue branch_queue_;
    size_t explored_node_count_;

    template<class OUT_STREAM>
    void print_status_bar_(OUT_STREAM& out)
    {
        if (verbosity == 0)
            return;
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
        for (size_t i = 0; i < callback_.num_separators(); ++i)
            std::cout << std::setw(20) << callback_.name(i) << " ";
        std::cout << "\n";
    }

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
        for (size_t i = 0; i < callback_.num_separators(); ++i)
        {

            std::cout 
                << std::setw(4) << callback_.num_calls(i) << " "
                << std::setw(7) << callback_.num_inequalities(i) << " " 
                << std::setw(7) << std::setprecision(5) << callback_.time(i) << " ";
        }
        if (verbosity >= 2)
            std::cout << "\n";
        else 
            std::cout << "    " << std::flush;

    }

    void add_constraint(const GRBTempConstr& constr)
    {
        current_node_->constraints_.emplace_back(current_node_->model.addConstr(constr));
        current_node_->num_non_basic_iterations_.push_back(0);
    }

    void add_constraint(const Inequality<int>& inequality)
    {
        GRBLinExpr expr;
        for (size_t i = 0; i < inequality.edges().size(); ++i)
        {
            size_t u = inequality.edges()[i][0];
            size_t v = inequality.edges()[i][1];
            expr += inequality.coefficients()[i] * vars_[u][v];
        }
        add_constraint(expr <= inequality.rhs());
    }

    size_t separate_()
    {   
        ++num_separation_calls_;

        VECTOR result(n_ * (n_-1) / 2);
        EPM result_map(n_, result);
        for (size_t i = 0; i < n_; ++i)
        for (size_t j = i+1; j < n_; ++j)
            result_map(i, j) = edge_value(i, j);

        return separate_continuous_(result_map);
    }

    size_t separate_continuous_(EPM edge_values)
    {
        auto inequalities = callback_(edge_values, stage_, current_node_->min_stage_);
        for (const auto& inequality : inequalities)
        {
            add_constraint(inequality);
        }
        return inequalities.size();
    }

    size_t remove_non_binding_constraints()
    {
        assert (current_node_->constraints_.size() + current_node_->depth == current_node_->model.get(GRB_IntAttr_NumConstrs));
        assert (current_node_->num_non_basic_iterations_.size() == current_node_->constraints_.size()); 
        size_t num_removed = 0;
        auto it_constr = current_node_->constraints_.cbegin();
        auto it_count = current_node_->num_non_basic_iterations_.begin();

        while ( it_constr != current_node_->constraints_.cend() ) 
        {   
            if (it_constr->get(GRB_IntAttr_CBasis) == 0) {
                ++(*it_count);
            }
            else {
                *it_count = 0;
            }
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

    void branch_()
    {
        if (num_integer() == n_ * (n_-1) / 2)
            throw std::runtime_error("current node is already integer!");

        // select a variable to branch on
        size_t i;
        size_t j; 
        select_branching_variable_(i, j);
        //
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

    void select_branching_variable_(size_t& i, size_t& j)
    {
        double largest_value = 0;
        size_t most_fractional_i;
        size_t most_fractional_j;
        double least_dist_to_one_half = 1;
        for (size_t a = 0; a < n_; ++a)
        for (size_t b = a+1; b < n_; ++b)
        {
            double x_ab = edge_value(a, b);
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

    void solve_current_node_lp_()
    {
        current_node_->model.update();
        tail_counter_ = 0;
        stage_ = 0;
        log_();
        while (true)
        {
            float_time_point lp_start = Time::now();
            double elapsed_time = (lp_start - optimization_start_time_).count();
            if (elapsed_time > time_limit_)
                break;
            current_node_->model.set(GRB_DoubleParam_TimeLimit, time_limit_ - elapsed_time);
            current_node_->model.optimize();
            float_time_point lp_end = Time::now();
            lp_time_ += (lp_end - lp_start).count();
            if (current_node_->model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
                throw std::runtime_error("Model is infeasible!");
            if (current_node_->model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                return;
            
            double new_objective = current_node_->model.get(GRB_DoubleAttr_ObjVal);
            double old_objective = current_node_->bound;
            current_node_->set_bound(new_objective);

            print_status_(std::cout);

            if (current_node_->model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
                break;

            if (stage_ > current_node_->min_stage_);
            {
                current_node_->min_stage_ = stage_;
            }

            // logging
            log_();

            // check if tailing of condition is met
            if (new_objective >= tail_threshold * old_objective)
            {
                ++tail_counter_;
                if (tail_counter_ >= max_tail_length)
                {
                    // If the minimum callback stage can be increase, increase the stage
                    if (current_node_->min_stage_ + 1 < callback_.num_stages())
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
                size_t num_removed = remove_non_binding_constraints();
            }

            if (activate_branching && current_node_->bound < best_objective_ + 1 - EPSILON)
            {
                break;
            }

            size_t n = separate_();
            if (n == 0)
            {
                break;
            }

            // if (current_node_->bound < root_node_.bound - EPSILON)
               //  return;
        }
        current_node_->lp_solved = true;
        ++explored_node_count_;
    }

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

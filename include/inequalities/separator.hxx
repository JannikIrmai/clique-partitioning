#include <inequalities/inequality.hxx>
#include <time.hxx>

#pragma once


namespace CP {

/**
 * Abstract class that represents a generic separation algorithm.
 * Separation algorithms are implemented as subclasses of this abstract class.
 * In particular each separation algorithm must implement the following functions:
 *  - name: returns the name of the separation algorithm
 *  - separate_: returns a list of violated inequalities
 * This base class support several helper functions for managing violated inequalities,
 * e.g. sorting and reducing the number of violated inequalities.
 * Moreover this class tracks some metrics for the separation algorithm like the 
 * number of times it was called, the number of inequalities that it found and the total
 * running time.
 */
template<class C, class EDGE_VALUE_MAP>
class AbstractSeparator
{
public:
    typedef Inequality<C> INEQUALITY;
    typedef std::vector<INEQUALITY> INEQUALITIES;

    AbstractSeparator(size_t n, size_t max_num = 0) : 
        n_(n),
        max_num_(max_num)
    {
        if (max_num_ == 0)
            max_num_ = n_*(n_-1)/2;
    }

    INEQUALITIES separate(const EDGE_VALUE_MAP& edge_values)
    {
        float_time_point start = Time::now();
        INEQUALITIES inequalities = separate_(edge_values);
        float_time_point end = Time::now();
        assert (inequalities.size() <= max_num_);
        total_num_ += inequalities.size();
        total_time_ += (end - start).count();
        ++num_calls_;
        return inequalities;
    }

    virtual std::string name() = 0;
    double total_time() {return total_time_;}
    size_t total_num() {return total_num_;}
    size_t num_calls() {return num_calls_;}

protected:
    size_t n_;  // number of nodes of the CP instance
    size_t max_num_;  // maximum number of inequalities per iteration
    size_t num_calls_ = 0;  // number of times the separator got called
    size_t total_num_ = 0;  // total number of inequalities over all iterations
    double total_time_ = 0;  // total separation time over all iterations

    double min_violation_ = 0.01;  // only inequalities whose violation is greater than this value are separated
    double min_violation_depth_ = 0.002;  // only inequalities whose violation depth is greater than this value are separated
    double min_relative_violation_depth_ = 0.5;  // only inequalities whose violation depth is greater than this value times the greatest violation depth of all violated inequalities are separated
    double max_parallel_ = 0.5;  // for two inequalities that are more parallel than this coefficient, only one is separated

    virtual INEQUALITIES separate_(const EDGE_VALUE_MAP& edge_values) = 0;

    // this method implements a heuristic of Sorensen (2020) for reducing the number of inequalities
    // based on their violation depth
    void sort_and_reduce_by_violation_depth(INEQUALITIES& inequalities) const
    {
        if (inequalities.size() == 0)
            return;
    
        // sort by violation depth
        std::stable_sort(inequalities.begin(), inequalities.end(),
            [&inequalities](const INEQUALITY& i1, const INEQUALITY& i2) {return i1.violation_depth() > i2.violation_depth();});
        // remove inequalities if their violation depth is less than min_violation_depth_ or if it is less than half of the largest violation depth
        double threshold = std::max(min_relative_violation_depth_*inequalities[0].violation_depth(), min_violation_depth_);
        for (size_t i = 0; i < inequalities.size(); ++i)
        {
            if (i >= max_num_ || inequalities[i].violation_depth() < threshold)
            {
                inequalities.resize(i);
                break;
            }
        }
    }

    // this method implements a heuristic of Sorensen (2020) for reducing the number of inequalities
    // based on how parallel they are
    void reduce_by_parallelism(INEQUALITIES& inequalities) const
    {
        std::vector<INEQUALITY> selected_inequalities;
        for (const INEQUALITY& ineq : inequalities)
        {
            bool skip = false;
            for (const INEQUALITY& ineq_selected : selected_inequalities)
            {
                double dot = ineq_selected.dot_product(ineq);
                double parallel_coeff = dot / (ineq_selected.euclidean() * ineq.euclidean());
                if (parallel_coeff > max_parallel_)
                {
                    skip =  true;
                    break;
                }
            }
            if (!skip)
            {
                selected_inequalities.push_back(ineq);
            }
        }
        inequalities = selected_inequalities;
    }

    // this method implements a heuristic for reducing the number of inequalities by
    // removing those inequalities that only contain non-zero coefficients for edges
    // that are already 'covered' by other inequalities.
    void sort_and_reduce_by_overlap(INEQUALITIES& inequalities) const
    {
        if (inequalities.size() == 0)
            return;
    
        // sort by violation depth
        std::stable_sort(inequalities.begin(), inequalities.end(),
            [&inequalities](const INEQUALITY& i1, const INEQUALITY& i2) {return i1.violation_depth() > i2.violation_depth();});

        std::vector<std::vector<size_t>> covered(n_, std::vector<size_t>(n_, 0));
        std::vector<INEQUALITY> selected_inequalities;
        for (const INEQUALITY& inequality : inequalities)
        {
            // check if there is an edge in the inequality that is not yet covered
            size_t is_covered = false;
            for (size_t i = 0; i < inequality.edges().size(); ++i)
            {
                size_t u = inequality.edges()[i][0];
                size_t v = inequality.edges()[i][1];
                if (inequality.coefficients()[i] < 0 && covered[u][v] == 1)
                {
                    is_covered = true;
                    break;
                }
                if (inequality.coefficients()[i] > 0 && covered[v][u] == 1)
                {
                    is_covered = true;
                    break;
                }
            }
            if (is_covered)
                continue;
            // mark all edge in the inequality as covered
            for (size_t i = 0; i < inequality.edges().size(); ++i)
            {
                size_t u = inequality.edges()[i][0];
                size_t v = inequality.edges()[i][1];
                if (inequality.coefficients()[i] < 0)
                    covered[u][v] = 1;
                else if (inequality.coefficients()[i] > 0)
                    covered[v][u] = 1;
                
            }
            selected_inequalities.push_back(inequality);
            if (selected_inequalities.size() >= max_num_)
                break;
        }
        inequalities = selected_inequalities;
    }
};


} // namespace CP

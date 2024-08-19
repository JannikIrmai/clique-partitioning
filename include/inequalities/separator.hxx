#include <inequalities/inequality.hxx>
#include <time.hxx>

#pragma once


namespace CP {


template<class C, class EDGE_VALUE_MAP>
class AbstractSeparator
{
public:
    typedef Inequality<C> INEQUALITY;
    typedef std::vector<INEQUALITY> INEQUALITIES;

    AbstractSeparator(size_t n, size_t max_num = 0) 
    : 
        n_(n),
        max_num_(max_num)
    {
        if (max_num_ == 0)
            max_num_ = n_*(n_-1)/2;
    }
    INEQUALITIES separate(EDGE_VALUE_MAP edge_values)
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
    double min_euclidean_violation_ = 0.002;  // only inequalities whose euclidean violation is greater than this value are separated
    double min_relative_euclidean_violation_ = 0.5;  // only inequalities whose euclidean violation is greater than this value times the greatest euclidean violation of all violated inequalities are separated
    double max_parallel_ = 0.5;  // for two inequalities that are more parallel than this coefficient, only one is separated

    virtual INEQUALITIES separate_(EDGE_VALUE_MAP edge_values) = 0;

    void sort_and_reduce_by_euclidean_violation(INEQUALITIES& inequalities) const
    {
        if (inequalities.size() == 0)
            return;
    
        // sort by euclidean violation
        std::stable_sort(inequalities.begin(), inequalities.end(),
            [&inequalities](const INEQUALITY& i1, const INEQUALITY& i2) {return i1.euclidean_violation() > i2.euclidean_violation();});
        // remove inequalities if their euclidean violation is less than 0.002 or if it is less than half of the largest euclidean violation
        double threshold = std::max(min_relative_euclidean_violation_*inequalities[0].euclidean_violation(), min_euclidean_violation_);
        for (size_t i = 0; i < inequalities.size(); ++i)
        {
            if (i >= max_num_ || inequalities[i].euclidean_violation() < threshold)
            {
                inequalities.resize(i);
                break;
            }
        }
    }

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

    void sort_and_reduce_by_overlap(INEQUALITIES& inequalities) const
    {
        if (inequalities.size() == 0)
            return;
    
        // sort by euclidean violation
        std::stable_sort(inequalities.begin(), inequalities.end(),
            [&inequalities](const INEQUALITY& i1, const INEQUALITY& i2) {return i1.euclidean_violation() > i2.euclidean_violation();});

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

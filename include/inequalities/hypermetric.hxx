#include <vector>
#include <map>
#include <cassert>
#include <numeric>
#include <inequalities/separator.hxx>
#include <inequalities/constants.hxx>
#include <inequalities/inequality.hxx>


namespace CP {


struct HyperMetric
{
    std::vector<int> b;
    int rhs;
    double violation;

    // values for computing euclidean violation
    int sum_of_squares;
    int sum_of_square_squares;
    double euclidean_violation;

    void print()
    {
        for (size_t i = 0; i < b.size(); ++i)
        if (b[i] != 0)
            std::cout << "b_" << i << " = " << b[i] << " ";
        std::cout << "\b, violation = " << violation << "\n";
    }

    template<class EDGE_VALUE_MAP>
    double compute_violation(EDGE_VALUE_MAP edge_values)
    {
        violation = - (double)rhs;
        for (size_t i = 0; i < b.size(); ++i)
        {
            if (b[i] == 0)
                continue;
            for (size_t j = i+1; j < b.size(); ++j)
            {
                if (b[j] == 0)
                    continue;
                violation -= edge_values(i, j) * b[i] * b[j];
            }
        }
        return violation;
    }

    int compute_rhs()
    {
        rhs = 0;
        for (auto bi : b)
        {
            rhs += bi * (bi - 1) / 2;
        }
        return rhs;
    }

    double compute_euclidean_violation()
    {
        sum_of_squares = 0;
        sum_of_square_squares = 0;
        for (auto bi : b)
        {
            int square = bi * bi;
            sum_of_squares += square;
            sum_of_square_squares += square * square;
        }
        double euclidean = std::sqrt((sum_of_squares * sum_of_squares - sum_of_square_squares) / 2);
        euclidean_violation = violation / euclidean;
        return euclidean_violation;
    }

    bool operator<(const HyperMetric& other) const
    {
        assert (b.size() == other.b.size());
        if (rhs < other.rhs)
            return true;
        else if (rhs > other.rhs)
            return false;
        if (sum_of_squares < other.sum_of_squares)
            return true;
        else if (sum_of_squares > other.sum_of_squares)
            return false;
        if (sum_of_square_squares < other.sum_of_square_squares)
            return true;
        else if (sum_of_square_squares > other.sum_of_square_squares)
            return false;
        for (size_t i = 0; i < b.size(); ++i)
        {
            if (b[i] < other.b[i])  
                return true;
            else if (b[i] > other.b[i])  
                return false;
        }
        return false;
    }
};


Inequality<int> hyper_metric_to_inequality(const HyperMetric& hyper_metric)
{
    // compute the number of non-zero coefficients in the inequality
    std::vector<std::array<size_t, 2>> edges;
    std::vector<int> coefficients;
    
    for (size_t i = 0; i < hyper_metric.b.size(); ++i)
    {
        {
            if (hyper_metric.b[i] == 0)
                continue;
            for (size_t j = i+1; j < hyper_metric.b.size(); ++j)
            {
                if (hyper_metric.b[j] == 0)
                    continue;
                edges.push_back({i, j});
                coefficients.push_back(- hyper_metric.b[i] * hyper_metric.b[j]);
            }
        }
    }

    return {edges, coefficients, hyper_metric.rhs, hyper_metric.violation};
}



template<class EDGE_VALUE_MAP>
class HyperMetricSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    HyperMetricSeparator(size_t n, size_t max_num = 0, size_t p = 10, bool do_swaps = true) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        delta_(n), best_delta_(n)
    {
        p_ = std::min(p, n);
        do_swaps_ = do_swaps;
    }

    std::string name()
    {
        return "HyperMetric";
    }

    std::vector<Inequality<int>> separate_(EDGE_VALUE_MAP edge_values)
    {
        std::vector<HyperMetric> hyper_metrics = separate_hyper_metric(edge_values);
        std::vector<Inequality<int>> inequalities(hyper_metrics.size());
        for (size_t i = 0; i < hyper_metrics.size(); ++i)
        {
            assert (std::abs(hyper_metrics[i].violation - hyper_metrics[i].compute_violation(edge_values)) < 1e-3);
            inequalities[i] = hyper_metric_to_inequality(hyper_metrics[i]);
            assert (std::abs(hyper_metrics[i].violation - (inequalities[i].evaluate(edge_values) - inequalities[i].rhs())) < 1e-3);
        }
        this->sort_and_reduce_by_euclidean_violation(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }

    std::vector<HyperMetric> separate_hyper_metric(EDGE_VALUE_MAP edge_values)
    {
        std::vector<HyperMetric> hyper_metrics;

        // iterate over all edges {i, j}
        for (size_t i = 0; i < this->n_; ++i)
        for (size_t j = i+1; j < this->n_; ++j)
        {
            // continue if the value of the edge {i, j} is integral
            double x_ij = edge_values(i, j);
            if (x_ij < 1 - this->min_violation_)
                add_hyper_metric(i, j, edge_values, hyper_metrics);
        }
        return hyper_metrics;
    }

    void add_hyper_metric(size_t u, size_t v, EDGE_VALUE_MAP edge_values, std::vector<HyperMetric>& hyper_metrics)
    {
        size_t n = edge_values.n();
        assert (n == this->n_);

        // initialize hyper metric inequality with b_u = 1, b_v = 1 and b_i = 0 for all u != i != v
        HyperMetric h;
        h.b = std::vector<int>(this->n_, 0);
        h.b[u] = 1;
        h.b[v] = 1;
        h.rhs = 0;
        h.sum_of_squares = 2;
        h.sum_of_square_squares = 2;
        h.violation = - edge_values(u, v);
        h.euclidean_violation = h.violation;

        size_t non_zeros = 2;
        HyperMetric best_h = h;
        size_t best_non_zeros = non_zeros;

        std::vector<int> change_direction(n, 0);
        
        // reset delta
        // delta[i] = sum_{j != i} b_j*x_{i, j}
        for (size_t i = 0; i < n; ++i)
        {
            delta_[i] = 0;
            if (i != u)
                delta_[i] += h.b[u] * edge_values(i, u);
            if (i != v)
                delta_[i] += h.b[v] * edge_values(i, v);
        }
        std::copy(delta_.begin(), delta_.end(), best_delta_.begin());

        // Compute the change in the violation when adding 
        // c_i to b_i for all i. It holds that
        //      - sum_{i < j} (b_i+c_i)*(b_j+c_j)*x_ij <= sum_i (b_i + c_i) * (b_i + c_i - 1) / 2
        // <=>  - sum_{i < j} b_i*b_j*x_ij - sum_i c_i delta_i - sum_{i < j} c_i*c_j*x_ij <= sum_i b_i*(b_i-1)/2 + \sum_i (b_i*c_i + c_i*(c_i-1)/2)
        // In particular, the left hand side of the hypermetric inequality changes by
        //      - sum_i c_i delta_i - sum_{i < j} c_i*c_j*x_ij 
        // and the right hand side changes by
        //      \sum_i (b_i*c_i + c_i*(c_i-1)/2)

        // compute change of the left hand side of the hypermetric inequality
        auto get_lhs_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            double lhs_change = 0;
            for (size_t i = 0; i < change.size(); ++i)
            {
                lhs_change -= change[i].second * delta_[change[i].first];
                for (size_t j = i+1; j < change.size(); ++j)
                {
                    lhs_change -= change[i].second * change[j].second * edge_values(change[i].first, change[j].first);
                }
            }
            return lhs_change;
        };

        // compute the change of the right hand side of the hypermetric inequality
        auto get_rhs_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            int rhs_change = 0;
            for (auto c : change)
            {
                rhs_change += h.b[c.first] * c.second + c.second * (c.second - 1) / 2;
            }
            return rhs_change;
        };

        // compute the change of the violation
        auto get_violation_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            return get_lhs_change(change) - get_rhs_change(change);
        };

        // The euclidean norm of the support vector of a hyper metric inequality is
        //      sqrt(sum_{i < j} (b_i*b_j)^2)
        // It holds that 
        //      2 * sum_{i < j} (b_i*b_j)^2 = (\sum_i b_i^2)^2 - \sum_i b_i^4 .
        auto get_squares_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            int new_sum_of_squares = h.sum_of_squares;
            int new_sum_of_square_squares = h.sum_of_square_squares;
            for (auto c : change)
            {
                new_sum_of_squares += c.second * (2 * h.b[c.first] + c.second);
                int b_i_plus_c_i = h.b[c.first] + c.second;
                int b_i_plus_c_i_squared = b_i_plus_c_i * b_i_plus_c_i;
                int b_i_squared = h.b[c.first] * h.b[c.first];
                new_sum_of_square_squares += b_i_plus_c_i_squared * b_i_plus_c_i_squared - b_i_squared * b_i_squared;
            }
            return std::pair<int, int>{new_sum_of_squares, new_sum_of_square_squares};
        };

        auto get_depth_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            double violation_change = get_violation_change(change);
            double new_violation = h.violation + violation_change;
            auto squares = get_squares_change(change);
            int new_squared_euclidean = (squares.first * squares.first - squares.second) / 2;
            if (new_squared_euclidean == 0)
                return -std::numeric_limits<double>::infinity();
            double new_euclidean = std::sqrt(new_squared_euclidean);
            double new_violation_depth = new_violation / new_euclidean;
            double depth_change = new_violation_depth - h.euclidean_violation;
            return depth_change;
        };

        auto make_change = [&] (const std::vector<std::pair<size_t, int>>& change)
        {
            double lhs_change = get_lhs_change(change);
            int rhs_change = get_rhs_change(change);
            double violation_change = lhs_change - rhs_change;
            h.violation += violation_change;
            h.rhs += rhs_change;
            auto squares = get_squares_change(change);
            h.sum_of_squares = squares.first;
            h.sum_of_square_squares = squares.second;
            int new_squared_euclidean = (squares.first * squares.first - squares.second) / 2;
            if (new_squared_euclidean == 0)
                throw std::runtime_error("Changing b leads to b == 0!");
            double new_euclidean = std::sqrt(new_squared_euclidean);
            h.euclidean_violation = h.violation / new_euclidean;
            for (auto c : change)
                h.b[c.first] += c.second;

            assert (h.rhs == h.compute_rhs());
            assert (std::abs(h.violation - h.compute_violation(edge_values)) < 1e-3 );
            assert (std::abs(h.euclidean_violation - h.compute_euclidean_violation()) < 1e-3 );

            // update delta
            for (size_t i = 0; i < n; ++i)
            {
                for (auto c : change)
                {
                    if (c.first != i)
                        delta_[i] += c.second * edge_values(i, c.first);
                }
            }
        };

        // compute the node i for which changing b[i] increases
        // the violation maximally
        auto get_best_node_change = [&] (size_t& best_i, int& best_change, double& best_violation_gain) 
        {
            best_violation_gain = -std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < n; ++i)
            {
                std::vector<int> changes;
                if (change_direction[i] == 0)
                    changes = {1, -1};
                else if (change_direction[i] > 0)
                    changes = {1};
                else 
                    changes = {-1};
            
                for (auto c : changes)
                {
                    double gain = get_violation_change({{i, c}});
                    if (gain > best_violation_gain)
                    {
                        best_i = i;
                        best_change = c;
                        best_violation_gain = gain;
                    }
                }
            }
        };
                        
        // ------- CONSTRUCTION PHASE ------------------------
        // increase/decrease the coefficients until there are p_
        // non-zero coefficients
        change_direction[u] = 1;
        change_direction[v] = 1;

        for (size_t k = 0; k < 2*p_; ++k)
        {
            size_t node;
            int change;
            double violation_gain;
            get_best_node_change(node, change, violation_gain);
            // if changing the value of that node would result in more than
            // p_ nodes with non-zero value, terminate the construction phase
            if (h.b[node] == 0)
                ++non_zeros;
            if (non_zeros > p_)
                break;
            
            make_change({{node, change}});
            change_direction[node] = change;
            // check if new hyper-metric is best
            if (h.euclidean_violation > best_h.euclidean_violation + EPSILON)
            {
                best_h = h;
                std::copy(delta_.begin(), delta_.end(), best_delta_.begin());
                best_non_zeros = non_zeros;
            }
        }

        assert (std::abs(best_h.violation - best_h.compute_violation(edge_values)) < 1e-3);
        assert (std::abs(best_h.euclidean_violation - best_h.compute_euclidean_violation()) < 1e-3);

        if (best_h.violation > this->min_violation_)
            hyper_metrics.push_back(best_h);
        double best_violation_depth_after_construction = best_h.euclidean_violation;
        
        if (best_h.violation <= -0.5)
            return;

        // --------- IMPROVEMENT PHASE -----------------------
        // try to improve the violation of the inequality by
        // increasing one node value while simultaneously
        // decreasing another node value.

        // Changing individual bi
        bool improvement = true;
        while (improvement)
        {
            improvement = false;
            h = best_h;
            std::copy(best_delta_.begin(), best_delta_.end(), delta_.begin());
            non_zeros = best_non_zeros;
            std::fill(change_direction.begin(), change_direction.end(), 0);

            for (size_t k = 0; k < p_; ++k)
            {
                double best_depth_gain = - std::numeric_limits<double>::infinity();
                size_t best_i;
                int best_change;
                int best_zero_change;
                for (size_t i = 0; i < n; ++i)
                {
                    std::vector<int> changes;
                    if (change_direction[i] < 0)
                        changes = {-1};
                    else if (change_direction[i] > 0)
                        changes = {1};
                    else
                        changes = {1, -1};
                    for (auto c : changes)
                    {
                        int zero_change = h.b[i] == 0 ? 1 : (h.b[i] == -c ? -1 : 0);
                        if (non_zeros + zero_change > p_)
                            continue;

                        double depth_gain = get_depth_change({{i, c}});
                        if (depth_gain > best_depth_gain)
                        {
                            best_depth_gain = depth_gain;
                            best_i = i;
                            best_zero_change = zero_change;
                            best_change = c;
                        }
                    }
                }
                make_change({{best_i, best_change}});
                change_direction[best_i] = best_change;
                non_zeros += best_zero_change;

                if (h.euclidean_violation > best_h.euclidean_violation + EPSILON)
                {
                    best_h = h;
                    std::copy(delta_.begin(), delta_.end(), best_delta_.begin());
                    best_non_zeros = non_zeros;
                    improvement = true;
                }
            }
        }

        if (!do_swaps_)
            return;

        h = best_h;
        std::copy(best_delta_.begin(), best_delta_.end(), delta_.begin());
        non_zeros = best_non_zeros;

        improvement = true;
        while (improvement)
        {
            improvement = false;
            for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
            {
                if (i == j)
                    continue;
                int i_zero_change = h.b[i] == 0 ? 1 : (h.b[i] == -1 ? -1 : 0);
                int j_zero_change = h.b[j] == 0 ? 1 : (h.b[j] == 1 ? -1 : 0);
                int zero_change = i_zero_change + j_zero_change;
                if (non_zeros + zero_change > p_)
                    continue;

                double swap_depth_gain = get_depth_change({{i, 1}, {j, -1}});
                if (swap_depth_gain > EPSILON)
                {
                    improvement = true;
                    make_change({{i, 1}, {j, -1}});
                    non_zeros += zero_change; 
                }
            }
        }
        
        if (h.euclidean_violation > best_violation_depth_after_construction && h.violation > this->min_violation_)
            hyper_metrics.push_back(h);
    }


private:

    size_t p_;
    bool do_swaps_;

    std::vector<double> delta_;
    std::vector<double> best_delta_;
};
    
} // namespace CP

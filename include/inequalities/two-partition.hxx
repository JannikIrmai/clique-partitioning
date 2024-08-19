#include <vector>
#include <cassert>
#include <numeric>
#include <random>

#include <inequalities/separator.hxx>
#include <inequalities/inequality.hxx>


namespace CP {


struct TwoPartition
{
    std::vector<size_t> s;
    std::vector<size_t> t;
    double violation;

    void print()
    {
        std::cout << "A = {";
        for (size_t a : s)
            std::cout << a << " ";
        std::cout << "\b}, B = {";
        for (size_t b : t)
            std::cout << b << " ";
        std::cout << "\b}, violation = " << violation << "\n";
    }

    template<class EDGE_VALUE_MAP>
    double compute_violation(EDGE_VALUE_MAP edge_values)
    {
        // compute violation of final two partition
        violation = - (double)std::min(s.size(), t.size());
        for (size_t i = 0; i < s.size(); ++i)
        for (size_t j = i+1; j < s.size(); ++j)
        {
            violation -= edge_values(s[i], s[j]);
        }
        for (size_t i = 0; i < t.size(); ++i)
        for (size_t j = i+1; j < t.size(); ++j)
        {
            violation -= edge_values(t[i], t[j]);
        }
        for (size_t i = 0; i < s.size(); ++i)
        for (size_t j = 0; j < t.size(); ++j)
        {
            violation += edge_values(s[i], t[j]);
        }
        return violation;
    }

};




Inequality<int> two_partition_to_inequality(const TwoPartition& two_partition)
{
    // compute the number of non-zero coefficients in the inequality
    size_t n = (two_partition.t.size() * (two_partition.t.size() - 1)) / 2 
        + (two_partition.s.size() * (two_partition.s.size() - 1)) / 2
        + two_partition.s.size() * two_partition.t.size();
    std::vector<std::array<size_t, 2>> edges(n);
    std::vector<int> coefficients(n);
    size_t idx = 0;

    // edges within s
    for (size_t i = 0; i < two_partition.s.size(); ++i)
    for (size_t j = i+1; j < two_partition.s.size(); ++j)
    {
        edges[idx] = {two_partition.s[i], two_partition.s[j]};
        coefficients[idx] = -1;
        ++idx;
    }

    // edges within t
    for (size_t i = 0; i < two_partition.t.size(); ++i)
    for (size_t j = i+1; j < two_partition.t.size(); ++j)
    {
        edges[idx] = {two_partition.t[i], two_partition.t[j]};
        coefficients[idx] = -1;
        ++idx;
    }

    // edges between s and t
    for (size_t i : two_partition.s)
    for (size_t j : two_partition.t)
    {
        edges[idx] = {i, j};
        coefficients[idx] = 1;
        ++idx;
    }

    assert (idx == n);

    return {edges, coefficients, (int)std::min(two_partition.s.size(), two_partition.t.size()), two_partition.violation};
}

/// separation heuristic by Sorensen (2020)

template<class EDGE_VALUE_MAP>
class TwoPartitionSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    TwoPartitionSeparator(size_t n, size_t max_num = 0, size_t p = 10) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num)
    {
        p_ = std::min(p, n);
    }

    std::string name()
    {
        return "TwoPartition";
    }

    std::vector<Inequality<int>> separate_(const EDGE_VALUE_MAP& edge_values)
    {
        size_t max_size = 0;        
        size_t total_size = 0;

        std::vector<TwoPartition> two_partitions = separate_two_partition(edge_values);
        std::vector<Inequality<int>> inequalities(two_partitions.size());
        for (size_t i = 0; i < two_partitions.size(); ++i)
        {
            size_t size = two_partitions[i].s.size() + two_partitions[i].t.size();
            total_size += size;
            if (size > max_size)
                max_size = size;
            inequalities[i] = two_partition_to_inequality(two_partitions[i]);
        }
        this->sort_and_reduce_by_euclidean_violation(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }


    std::vector<TwoPartition> separate_two_partition(const EDGE_VALUE_MAP& edge_values)
    {
        std::vector<TwoPartition> two_partitions;

        // iterate over all edges {i, j}
        for (size_t i = 0; i < this->n_; ++i)
        for (size_t j = i+1; j < this->n_; ++j)
        {
            // continue if the value of the edge {i, j} is integral
            double x_ij = edge_values(i, j);
            if (x_ij <= this->min_violation_ || (1 - x_ij) <= this->min_violation_)
                continue;

            add_two_partition(i, j, edge_values, two_partitions);
        }

        return two_partitions;
    }

    void add_two_partition(size_t a, size_t b, const EDGE_VALUE_MAP& edge_values, std::vector<TwoPartition>& two_partitions)
    {
        // ------- CONSTRUCTION PHASE ------------------------
        // start with the partition A = {a} and B = {b} and
        // add nodes to either A or B until |A \cup B| = p_.
        // For each added node k, we track how much the violation
        // of the inequality changes by adding that node.
        std::vector<double> violation_gains = {-1, edge_values(a, b)};
        std::vector<size_t> nodes_in_partition = {a, b};
        std::vector<size_t> A_or_B = {1, 0};  // characteristic vector indicating,
        // whether a node in nodes_in_partition is in A or B.
        size_t size_A = 1;
        size_t size_B = 1;

        double inf = std::numeric_limits<double>::infinity();
        // delta_A[i] is the sum of edge value x_{a,i} for 
        // all a in A. delta_B is analogue
        std::vector<double> delta_A(this->n_);
        std::vector<double> delta_B(this->n_);
        for (size_t i = 0; i < this->n_; ++i)
        {
            if (i == a || i == b) {
                delta_A[i] = -inf;
                delta_B[i] = -inf;
            } else {
                delta_A[i] = edge_values(a, i);
                delta_B[i] = edge_values(b, i);
            }
        }

        for (size_t k = 2; k < p_; ++k)
        {
            double max_gain = -inf;
            size_t max_gain_node;
            bool max_gain_for_A;
            for (size_t i = 0; i < this->n_; ++i)
            {
                if (delta_A[i] == -inf)
                    continue;  // i is in A or B
                double gain_A = delta_B[i] - delta_A[i];
                if (size_A < size_B)
                    gain_A -= 1;
                if (gain_A > max_gain)
                {
                    max_gain = gain_A;
                    max_gain_node = i;
                    max_gain_for_A = true;
                }
                double gain_B = delta_A[i] - delta_B[i];
                if (size_B < size_A)
                    gain_B -= 1;
                if (gain_B > max_gain)
                {
                    max_gain = gain_B;
                    max_gain_node = i;
                    max_gain_for_A = false;
                }
            }
            // add max_gain_node to the respective set
            violation_gains.push_back(max_gain);
            nodes_in_partition.push_back(max_gain_node);
            A_or_B.push_back(max_gain_for_A);
            if (max_gain_for_A)
                ++size_A;
            else
                ++size_B;
            // std::cout << "adding " << max_gain_node << " to " << (max_gain_for_A? "A" : "B") << "\n";

            // update deltas
            for (size_t i = 0; i < this->n_; ++i)
            {
                if (i == max_gain_node)
                {
                    delta_A[i] = -inf;
                    delta_B[i] = -inf;
                    continue;
                }
                if (max_gain_for_A)
                    delta_A[i] += edge_values(max_gain_node, i);
                else
                    delta_B[i] += edge_values(max_gain_node, i);
            }
        }

        // get best partition along sequence
        double violation = 0;
        double max_violation = 0;
        double max_euclidean_violation = -inf;
        size_t max_idx;
        for (size_t k = 0; k < p_; ++k)
        {
            violation += violation_gains[k];
            double euclidean_violation = violation / std::sqrt((k * (k+1)) / 2);
            if (euclidean_violation > max_euclidean_violation)
            {
                max_violation = violation;
                max_euclidean_violation = euclidean_violation;
                max_idx = k;
            }
        }
        // extract the partitions up to max_idx
        std::vector<size_t> A;
        std::vector<size_t> B;
        for (size_t k = 0; k <= max_idx; ++k)
        {
            if (A_or_B[k] == 1)
                A.push_back(nodes_in_partition[k]);
            else
                B.push_back(nodes_in_partition[k]);
        }

        
        // add two partition if it is violated
        if (max_violation > this->min_violation_)
        {
            TwoPartition tp = {A, B, max_violation};
            assert (std::abs(max_violation - tp.compute_violation(edge_values)) < 1e-3);
            two_partitions.push_back({A, B, max_violation});
        }

        // ------- Improvement phase ---------------
        
        // try removing nodes from A and B
        while (A.size() > 1 && B.size() > 1)
        {
            bool removed_node = false;
            // try to remove nodes from A
            for (auto ait = A.begin(); ait != A.end(); ++ait)
            {
                double gain = 0;
                for (size_t b : B)
                    gain -= edge_values(*ait, b);
                for (size_t a2 : A)
                    if (a2 != *ait)
                        gain += edge_values(*ait, a2);
                if (A.size() <= B.size())
                    gain += 1;
                if (gain > 0);
                {
                    // std::cout << "removing " << *ait << " from A\n";
                    std::swap(*ait, A.back());
                    A.pop_back();
                    removed_node = true;
                    break;
                }
            }

            if (removed_node)
                continue;

            // try to remove nodes from B
            for (auto bit = B.begin(); bit != B.end(); ++bit)
            {
                double gain = 0;
                for (size_t a : A)
                    gain -= edge_values(*bit, a);
                for (size_t b2 : B)
                    if (b2 != *bit)
                        gain += edge_values(*bit, b2);
                if (B.size() <= A.size())
                    gain += 1;
                if (gain > 0);
                {
                    // std::cout << "removing " << *bit << " from B\n";
                    std::swap(*bit, B.back());
                    B.pop_back();
                    removed_node = true;
                    break;
                }
            }
            
            if (! removed_node)
                break;
        }
        
        // try to swap nodes in A and B with nodes not in A nor in B
        bool swap = true;
        while (swap)
        {
            swap = false;
            for (size_t i = 0; i < this->n_; ++i)
            {
                // continue if i is in A or B
                if (std::find(A.begin(), A.end(), i) != A.end())
                    continue;
                if (std::find(B.begin(), B.end(), i) != B.end())
                    continue;
            
                // try swapping i with a node in A
                for (auto ait = A.begin(); ait != A.end(); ++ait)
                {
                    double gain = 0;
                    for (size_t b : B)
                    {
                        gain -= edge_values(*ait, b);
                        gain += edge_values(i, b);
                    }
                    for (size_t a2 : A)
                    {
                        if (a2 != *ait)
                        {
                            gain += edge_values(*ait, a2);
                            gain -= edge_values(i, a2);
                        }
                    }
                    if (gain > 1e-6)
                    {
                        *ait = i;
                        swap = true;
                        break;
                    }
                }
                if (swap)
                    break;

                // try swapping i with a node in B
                for (auto bit = B.begin(); bit != B.end(); ++bit)
                {
                    double gain = 0;
                    for (size_t a : A)
                    {
                        gain -= edge_values(*bit, a);
                        gain += edge_values(i, a);
                    }
                    for (size_t b2 : B)
                    {
                        if (b2 != *bit)
                        {
                            gain += edge_values(*bit, b2);
                            gain -= edge_values(i, b2);
                        }
                    }
                    if (gain > 1e-6)
                    {
                        *bit = i;
                        swap = true;
                        break;
                    }
                }
                if (swap)
                    break;
            }            
        }

        // compute violation of final two partition
        TwoPartition tp = {A, B};
        violation = tp.compute_violation(edge_values);

        if (violation > max_violation && violation > this->min_violation_)
        {
            two_partitions.push_back(tp);
        }
    }


private:

    size_t p_;
};

} // namespace CP
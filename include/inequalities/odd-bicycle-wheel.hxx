#pragma once

#include <vector>
#include <array>
#include <limits>
#include <queue>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <inequalities/inequality.hxx>
#include <inequalities/separator.hxx>


namespace CP {


struct BicycleWheel
{
    std::vector<size_t> wheel;
    std::array<size_t, 2> center;
    double violation;
};


template<class EDGE_VALUE_MAP>
class OddBicycleWheelSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    OddBicycleWheelSeparator(size_t n, size_t max_num = 0) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        predecessor(n),
        distance(n)
    {}
    
    std::string name() { return "OddBicycleWheel"; }

    std::vector<Inequality<int>> separate_(const EDGE_VALUE_MAP& edge_values)
    {
        std::vector<BicycleWheel> wheels = separate_bicycle_wheel(edge_values);
        std::vector<Inequality<int>> inequalities(wheels.size());
        for (size_t i = 0; i < wheels.size(); ++i)
        {
            size_t num_vars = 3*wheels[i].wheel.size() + 1;
            std::vector<std::array<size_t,2>> edges(num_vars);
            std::vector<int> coefficients(num_vars);
            edges[0] = wheels[i].center;
            coefficients[0] = -1;
            for (size_t j = 0; j < wheels[i].wheel.size(); ++j)
            {
                edges[3*j+1] = {wheels[i].center[0], wheels[i].wheel[j]};
                coefficients[3*j+1] = 1;
                edges[3*j+2] = {wheels[i].center[1], wheels[i].wheel[j]};
                coefficients[3*j+2] = 1;
                edges[3*j+3] = {wheels[i].wheel[j], wheels[i].wheel[(j+1) % wheels[i].wheel.size()]};
                coefficients[3*j+3] = -1;
            }
            inequalities[i] = {edges, coefficients, (int)wheels[i].wheel.size() - 1, wheels[i].violation};
        }
        this->sort_and_reduce_by_euclidean_violation(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }

    /**
     * @brief 
        Select two center nodes i, j.
        Consider a bipartite graph that contains two copies of 
        each node of the original graph except i and j.
        We denote the nodes of the bipartite graph by (k, 0)
        and (k, 1) for all k in {0,...,n-1} \ {i, j}.
        For every edge {k1, k2} with k1, k2 in {0,...,n-1} \ {i, j}, 
        k1 != k2 of the original graph, the bipartite graph contains 
        the edges {(k1, 0), (k2, 1)} and {(k2, 0), (k2, 1)} with 
        weights 1 + x_k1k2 - 1/2 (x_ik1 + x_ik2 + x_jk1 + x_jk2).
        Then there exists a violated odd wheel inequality with 
        center i and j iff there exists a (k1, 0) - (k1, 1) - path of 
        total weight strict less than (1 - x_ij) in the bipartite graph.
    * 
    * @tparam EDGE_PROPERTY_MAP 
    * @param edge_values 
    * @return std::vector<BicycleWheel> 
    */
    std::vector<BicycleWheel> separate_bicycle_wheel(const EDGE_VALUE_MAP& edge_values)
    {
        size_t n = edge_values.n();
        assert (n == this->n_);

        std::vector<BicycleWheel> wheels;
        // iterate over all potential centers
        for (size_t i = 0; i < n; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            if (edge_values(i, j) >= 1 - this->min_violation_)
                continue;  // inequality cannot be violated

            auto weight = [&](size_t s, size_t t) {
                return std::max(0.0, 1.0 + edge_values(s, t) - (edge_values(i, s) + edge_values(i, t) + edge_values(j, s) + edge_values(j, t)) / 2);
            };

            auto is_zero_spoke = [&](size_t s) {
                return (edge_values(i, s) <= this->min_violation_ || 
                        edge_values(j, s) <= this->min_violation_);
            };

            auto is_one_spoke = [&](size_t s) {
                return (edge_values(i, s) >= 1 - this->min_violation_ || 
                        edge_values(j, s) >= 1 - this->min_violation_);
            };
            
            for (size_t k = 0; k < n; ++k)
            {
                if (i == k || j == k)
                    continue;

                if (is_one_spoke(k))
                    continue; // spokes with value 1 cannot violate inequality

                // search for a shortest paths from k-even
                std::fill(distance.begin(), distance.end(), std::array<double, 2>({inf, inf}));
                
                distance[k][0] = 0;
                std::priority_queue<Vertex> queue;
                queue.push({k, 0, distance[k][0]});

                while (!queue.empty())
                {
                    Vertex v = queue.top();
                    queue.pop();
                    if (distance[v.id][v.even] < v.distance) 
                        continue;  // v is a copy with outdated distance
                    // iterate over all neighbors
                    for (size_t t = k; t < n; ++t)
                    {
                        if (t == i || t == j || t == v.id)
                            continue;
                        if (is_one_spoke(t))
                            continue; // spokes with value 1 cannot violate inequality
                        if (is_zero_spoke(t) && (is_zero_spoke(v.id) || is_one_spoke(v.id)))
                            continue; // two consecutive integer spokes cannot violate inequality

                        double alt_distance = distance[v.id][v.even] + weight(v.id, t);
                        size_t t_even = 1-v.even;  // if v was odd, t will be even and vice versa
                        if (distance[t][t_even] <= alt_distance)
                            continue;  // t was already visited and has less distance
                        if (alt_distance >= 1 - edge_values(i, j) - this->min_violation_)
                            continue;
                        distance[t][t_even] = alt_distance;
                        predecessor[t][t_even] = v.id;
                        queue.push({t, t_even, alt_distance});
                    }
                }

                double violation = (1.0 - edge_values(i, j)) - distance[k][1];
                if (violation <= this->min_violation_)
                    continue;
                                
                BicycleWheel wheel;
                wheel.center = {i, j};
                wheel.violation = violation;
                size_t end = k;
                size_t end_even = 1;
                // reconstruct the path
                wheel.wheel = {k};
                while (predecessor[end][end_even] != k || end_even == 0){
                    end = predecessor[end][end_even];
                    end_even = 1 - end_even;
                    wheel.wheel.push_back(end);
                }
                // assert that the wheel is indeed odd.
                assert (wheel.wheel.size() % 2 == 1);

                // mark the reverse wheel as found
                wheels.push_back(wheel);
            }
        }
        return wheels;
    }


private:
    std::vector<std::array<size_t, 2>> predecessor;
    std::vector<std::array<double, 2>> distance;
    double inf = std::numeric_limits<double>::infinity();

    struct Vertex
    {
        size_t id;
        size_t even;
        double distance;
        bool operator<(const Vertex& other) const { 
            return distance > other.distance; 
        }
    };
};

} // namespace CP

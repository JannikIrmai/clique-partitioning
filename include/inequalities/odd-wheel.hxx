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


struct Wheel
{
    std::vector<size_t> wheel;
    size_t center;
    double violation;
};


template<class EDGE_VALUE_MAP>
class OddWheelSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    OddWheelSeparator(size_t n, size_t max_num = 0) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        predecessor(n),
        distance(n),
        reverse_found(n)
    {}
    
    std::string name() { return "OddWheel"; }

    std::vector<Inequality<int>> separate_(const EDGE_VALUE_MAP& edge_values)
    {
        std::vector<Wheel> wheels = separate_wheel(edge_values);
        std::vector<Inequality<int>> inequalities(wheels.size());
        for (size_t i = 0; i < wheels.size(); ++i)
        {
            std::vector<std::array<size_t,2>> edges(2*wheels[i].wheel.size());
            std::vector<int> coefficients(2*wheels[i].wheel.size());
            for (size_t j = 0; j < wheels[i].wheel.size(); ++j)
            {
                edges[2*j] = {wheels[i].center, wheels[i].wheel[j]};
                coefficients[2*j] = 1;
                edges[2*j+1] = {wheels[i].wheel[j], wheels[i].wheel[(j+1) % wheels[i].wheel.size()]};
                coefficients[2*j+1] = -1;
            }
            inequalities[i] = {edges, coefficients, (int)wheels[i].wheel.size() / 2, wheels[i].violation};
        }
        this->sort_and_reduce_by_violation_depth(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }

    /**
     * @brief 
        Select a center node i.
        Consider a bipartite graph that contains two copies of 
        each node of the original graph except i.
        We denote the nodes of the bipartite graph by (j, 0)
        and (j, 1) for all j = 0,...,n-1 with j != i.
        For every edge {j, k} with j, k in {0,...,n} \ {i}, 
        j != k of the original graph, the bipartite graph contains 
        the edges {(j, 0), (k, 1)} and {(k, 0), (j, 1)} with 
        weights 1/2 + x_jk - 1/2 (x_ij + x_ik).
        Then there exists a violated odd wheel inequality with 
        center i iff there exists a (j, 0) - (j, 1) - path of 
        total weight strict less than 1/2 in the bipartite graph.
        NOTE: shortest paths in the bipartite graph only correspond
        to odd walks in the original graph (i.e. the same edge may
        be used twice, e.g. the path 1, 2, 3, 4, 2, 1). But this
        is not a problem as odd wheels are also valid in that case.
    * 
    * @tparam EDGE_PROPERTY_MAP 
    * @param edge_values 
    * @return std::vector<Wheel> 
    */
    std::vector<Wheel> separate_wheel(const EDGE_VALUE_MAP& edge_values)
    {
        size_t n = edge_values.n();
        assert (n == this->n_);

        std::vector<Wheel> wheels;
        // iterate over all potential centers
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                if (i == j)
                    continue;

                if (edge_values(i, j) <= this->min_violation_ || edge_values(i, j) >= 1 - this->min_violation_)
                    continue;  // odd wheel inequality cannot be violated for integer spokes

                // search for a shortest paths from j-even
                auto weight = [&edge_values, &i](size_t s, size_t t) {
                    return std::max(0.0, 0.5 + edge_values(s, t) - (edge_values(i, s) + edge_values(i, t)) / 2);
                };

                std::fill(distance.begin(), distance.end(), std::array<double, 2>({inf, inf}));
                std::fill(reverse_found.begin(), reverse_found.end(), false);
                
                distance[j][0] = 0;
                std::priority_queue<Vertex> queue;
                queue.push({j, 0, distance[j][0]});

                while (!queue.empty())
                {
                    Vertex v = queue.top();
                    queue.pop();
                    if (distance[v.id][v.even] < v.distance)
                        continue;  // v is a copy with outdated distance
                    // iterate over all neighbors
                    for (size_t t = j + 1; t < n; ++t)
                    {
                        if (t == i || t == v.id)
                            continue;  // i is center and s is current node
                        if (edge_values(i, t) <= this->min_violation_ || edge_values(i, t) >= 1 - this->min_violation_)
                            continue;  // odd wheel inequality cannot be violated for integer spokes
                        double alt_distance = distance[v.id][v.even] + weight(v.id, t);
                        size_t t_even = 1-v.even;  // if v was odd, t will be even and vice versa
                        if (distance[t][t_even] <= alt_distance)
                            continue;  // t was already visited and has less distance
                        if (alt_distance >= 0.5 - this->min_violation_)
                            continue;
                        distance[t][t_even] = alt_distance;
                        predecessor[t][t_even] = v.id;
                        queue.push({t, t_even, alt_distance});
                    }
                }

                for (size_t k = j+1; k < n; ++k)
                {
                    if (i == k)
                        continue;
                    if (reverse_found[k])
                        continue;
                    // the right hand side of the inequality is
                    // the weight of the edge jk in the auxiliary graph
                    // plus the lengths of the shortest path from j even to k even
                    double rhs = weight(j, k) + distance[k][0];

                    double violation = 0.5 - rhs;

                    if (violation < this->min_violation_)
                        continue;
                                    
                    Wheel wheel;
                    wheel.center = i;
                    wheel.violation = violation;
                    size_t end = k;
                    size_t end_even = 0;
                    // reconstruct the path
                    wheel.wheel = {j, k};
                    while (predecessor[end][end_even] != j || end_even == 0){
                        end = predecessor[end][end_even];
                        end_even = 1 - end_even;
                        wheel.wheel.push_back(end);
                    }
                    // assert that the wheel is indeed odd.
                    assert (wheel.wheel.size() % 2 == 1);

                    // mark the reverse wheel as found
                    reverse_found[end] = true;
                    wheels.push_back(wheel);
                }
            }
        }
        return wheels;
    }


private:
    std::vector<std::array<size_t, 2>> predecessor;
    std::vector<std::array<double, 2>> distance;
    std::vector<bool> reverse_found;
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

#include <vector>
#include <array>
#include <limits>
#include <queue>
#include <algorithm>
#include <numeric>

#include <inequalities/constants.hxx>
#include <inequalities/inequality.hxx>
#include <inequalities/separator.hxx>

namespace CP {


struct Cycle
{
    std::vector<size_t> cycle;
    double violation;
};


template<class EDGE_VALUE_MAP>
class HalfChordedOddCycleSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    HalfChordedOddCycleSeparator(size_t n, size_t max_num = 0) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        predecessor(n, std::vector<std::array<size_t, 2>>(n)),
        distance(n, std::vector<std::array<double, 2>>(n))
    {}
    
    std::string name() { return "HCO-Cycle"; }

    std::vector<Inequality<int>> separate_(EDGE_VALUE_MAP edge_values)
    {
        std::vector<Cycle> cycles = separate_cycles(edge_values);
        std::vector<Inequality<int>> inequalities(cycles.size());
        for (size_t i = 0; i < cycles.size(); ++i)
        {
            std::vector<std::array<size_t,2>> edges(2*cycles[i].cycle.size());
            std::vector<int> coefficients(2*cycles[i].cycle.size());
            for (size_t j = 0; j < cycles[i].cycle.size(); ++j)
            {
                edges[2*j] = {cycles[i].cycle[j], cycles[i].cycle[(j+1) % cycles[i].cycle.size()]};
                coefficients[2*j] = -1;
                edges[2*j+1] = {cycles[i].cycle[j], cycles[i].cycle[(j+2) % cycles[i].cycle.size()]};
                coefficients[2*j+1] = 1;
            }
            inequalities[i] = {edges, coefficients, (int)cycles[i].cycle.size() - 3, cycles[i].violation};
        }
        this->sort_and_reduce_by_euclidean_violation(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }

    std::vector<Cycle> separate_cycles(EDGE_VALUE_MAP edge_values)
    {
        size_t n = edge_values.n();
        assert (n == this->n_);

        // search for a shortest paths from j-even
        auto weight = [&edge_values](size_t i, size_t j, size_t k) {
            return std::max(0.0, edge_values(i, j) + (1 - edge_values(i, k)));
        };

        std::vector<Cycle> cycles;

        // iterate over all pairs of nodes
        for (size_t i = 0; i < n; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            if (edge_values(i, j) >= 1 - this->min_violation_)
                continue;  // cannot be violated

            // search for a shortest paths from the pair ij-even to ij-odd in the auxiliary graph

            // set the distances to all pairs to infinity
            std::fill(distance.begin(), distance.end(), std::vector<std::array<double,2>>(n, {inf, inf}));

            float_time_point start = Time::now();
            distance[i][j][0] = 0;
            std::priority_queue<Edge> queue;
            queue.push({i, j, 0, 0});

            while (!queue.empty())
            {   
                Edge e = queue.top();
                queue.pop();
                if (distance[e.i][e.j][e.even] < e.distance)
                    continue;  // e is a copy with outdated distance
                if (e.i == i && e.j == j && e.even == 1)
                    break;  // found the shortest path
                // iterate over all neighbors
                for (size_t k = 0; k < n; ++k)
                {
                    if (k == e.i || k == e.j)
                        continue; 
                    if (edge_values(e.i, k) <= EPSILON)
                        continue; // cannot be violated
                    double alt_distance = e.distance + weight(e.i, e.j, k);
                    size_t jk_even = 1-e.even;  // if ij was odd, jk will be even and vice versa
                    if (distance[e.j][k][jk_even] <= alt_distance + EPSILON)
                        continue;  // jk was already visited and has less distance
                    if (alt_distance >= 3 - this->min_violation_)
                        continue; // path is to long
                    distance[e.j][k][jk_even] = alt_distance;
                    predecessor[e.j][k][jk_even] = e.i;
                    queue.push({e.j, k, jk_even, alt_distance});
                }
            }

            // check if there exists a cycle of length less than 3
            double violation = 3 - distance[i][j][1];
            if (violation < this->min_violation_)
                continue;
                
            // extract shortest cycle
            Cycle cycle;
            cycle.violation = violation;
            cycle.cycle = {};
            size_t y = j;
            size_t x = i;
            size_t even = 1;
            while (even != 0 || x != i || y != j)
            {
                cycle.cycle.push_back(y);
                size_t new_x = predecessor[x][y][even];
                even = 1 - even;
                y = x;
                x = new_x;
            }
            // assert that the cycle has odd length
            assert (cycle.cycle.size() % 2 == 1);

            cycles.push_back(cycle);
        }
        return cycles;
    }


private:
    std::vector<std::vector<std::array<size_t, 2>>> predecessor;
    std::vector<std::vector<std::array<double, 2>>> distance;
    double inf = std::numeric_limits<double>::infinity();

    struct Edge
    {
        size_t i;
        size_t j;
        size_t even;
        double distance;
        bool operator<(const Edge& other) const { 
            return distance > other.distance; 
        }
    };
};

} // namespace CP

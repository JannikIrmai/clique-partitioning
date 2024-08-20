#include <vector>
#include <array>
#include <limits>
#include <queue>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <tuple>

#include <inequalities/inequality.hxx>
#include <inequalities/separator.hxx>


namespace CP {

/**
 * Class for separating {0,1/2}-Chvatal-Gomory cuts as established by Caprara and Fischetti (1996)
 */
template<class EDGE_VALUE_MAP>
class ChvatalGomorySeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    ChvatalGomorySeparator(size_t n, size_t max_num = 0) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        predecessor(n, std::vector<std::array<size_t, 2>>(n)),
        distance(n)
    {}
    
    std::string name() { return "ChvatalGomory"; }

    std::vector<Inequality<int>> separate_(const EDGE_VALUE_MAP& edge_values)
    {
        size_t n = edge_values.n();
        assert (n == this->n_);

        // The weight of the edges {ij, jk} in the auxiliary graph are the half of the slack of the relaxed
        // triangle inequality 0.5 - 0.5*(x_ij + x_jk) + x_ik. Clearly these weights are non-negative
        // if the triangle inequalities are satisfied.
        auto weight = [&edge_values](size_t i, size_t j, size_t k) {
            return std::max(0.0, 0.5 - 0.5*(edge_values(i, j) + edge_values(j, k)) + edge_values(i, k));
        };
        // there exists a violated inequality iff there exists a cycle of odd length in
        // the auxiliary graph with weight strictly less than 0.5.

        std::vector<Inequality<int>> inequalities;
        for (size_t i = 0; i < n; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            if (edge_values(i, j) >= 1 - this->min_violation_)
                continue;  // inequality cannot be violated

            // set all distances to infinity
            std::fill(distance.begin(), distance.end(), std::vector<double>(n, inf));

            // search for shortest paths starting from ij
            distance[i][j] = 0.0;
            std::priority_queue<Vertex> queue;
            queue.push({i, j, 0.0});

            while (!queue.empty())
            {
                Vertex v = queue.top();
                queue.pop();
                if (distance[v.i][v.j] < v.distance)
                    continue;  // v is a copy with outdated distance
                // iterate over all neighbors k >= i because for smaller
                // k all edges were checked in a previous iteration
                for (size_t k = i; k < n; ++k)
                {
                    if (k == v.i || k == v.j)
                        continue;
                    if (edge_values(v.i, k) >= 1 - this->min_violation_ ||
                        edge_values(v.j, k) >= 1 - this->min_violation_)
                        continue;  // inequality cannot be violated
                    for (size_t a = 0; a < 2; ++a)
                    {
                        size_t vi, vj;
                        if (a == 0)
                            vi = v.i, vj = v.j;
                        else
                            vi = v.j, vj = v.i;
                        if ((vi == i && k < j) || (k == i && vi < j))
                            continue;  // these edges were checkt in a previous iteration
                        // path to vi,k
                        double alt_distance = v.distance + weight(vj, vi, k);
                        size_t s, t;
                        if (v.i < v.j) {  // was even
                            if (vi < k)
                                s = k, t = vi;
                            else 
                                s = vi, t = k;
                        } else {  // was odd
                            if (vi < k)
                                s = vi, t = k;
                            else 
                                s = k, t = vi;
                        }            
                        if (distance[s][t] <= alt_distance)
                            continue;  // already visited and has less distance
                        if (alt_distance >= 0.5 - this->min_violation_)
                            continue;  // distance to long, inequality not violated
                        distance[s][t] = alt_distance;
                        predecessor[s][t] = {v.i, v.j};
                        queue.push({s, t, alt_distance});
                    }
                }
            }

            // extract shortest cycle from ij to ji
            double violation = 0.5 - distance[j][i];
            if (violation <= this->min_violation_)
                continue;
            std::array<size_t, 2> end = {j, i};
            std::vector<std::array<size_t,2>> edges;
            std::vector<int> coefficients;
            while (end[0] != i || end[1] != j)
            {
                edges.push_back(end);
                coefficients.push_back(1);
                std::array<size_t, 2> new_end = predecessor[end[0]][end[1]];
                if (end[0] == new_end[0])
                    edges.push_back({end[1], new_end[1]});
                else if (end[0] == new_end[1])
                    edges.push_back({end[1], new_end[0]});
                else if (end[1] == new_end[0])
                    edges.push_back({end[0], new_end[1]});
                else if (end[1] == new_end[1])
                    edges.push_back({end[0], new_end[0]});
                else
                    throw std::runtime_error("This cannot happen!");
                coefficients.push_back(-1);
                end = new_end;
            }
            assert (edges.size() % 2 == 0);  // assert that cycle is odd (i.e. even number of edges)
            inequalities.push_back({edges, coefficients, (int)edges.size() / 4, violation});
        }

        this->sort_and_reduce_by_violation_depth(inequalities);
        this->reduce_by_parallelism(inequalities);
        return inequalities;
    }

private:
    std::vector<std::vector<std::array<size_t, 2>>> predecessor;
    std::vector<std::vector<double>> distance;
    double inf = std::numeric_limits<double>::infinity();

    struct Vertex
    {
        size_t i;
        size_t j;
        // if i < j then the vertex is even, else odd
        double distance;
        bool operator<(const Vertex& other) const { 
            return distance > other.distance; 
        }
    };
};

} // namespace CP

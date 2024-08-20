#include <vector>
#include <array>
#include <set>
#include <numeric>
#include <algorithm>

#include <inequalities/separator.hxx>
#include <andres/partition.hxx>

#pragma once


namespace CP {

/**
 * This class implements an algorithm for separating triangle inequalities
 */
template<class EDGE_VALUE_MAP>
class TriangleSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP>
{
public :
    TriangleSeparator(size_t n, size_t max_num = 0) : AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num) {}

    std::string name() { return "Triangle"; }

    /**  
     * Compute violated triangle inequalities.
     * For each edge at most one triangle containing that edge is returned.
     */
    std::vector<Inequality<int>> separate_(const EDGE_VALUE_MAP& edge_values)
    {
        size_t n = edge_values.n();
        assert (edge_values.n() == this->n_);

        struct Triangle
        {
            std::array<size_t, 3> triangle;
            double violation = 0;
            bool operator<(const Triangle& other) const
            {
                if (triangle[0] == other.triangle[0]) {
                    if (triangle[1] == other.triangle[1])
                        return triangle[2] < other.triangle[2];
                    else
                        return triangle[1] < other.triangle[1];
                } else {
                    return triangle[0] < other.triangle[0];
                }
            }
        };

        std::vector<std::vector<Triangle>> edge2triangle(n, std::vector<Triangle>(n));
        std::vector<size_t> neighbors(n);
        std::iota(neighbors.begin(), neighbors.end(), 0);
        // iterate over all nodes i and search for violated triangle
        // inequalities of the form x_ij + x_ik - x_jk <= 1
        for (size_t i = 0; i < n; ++i)
        {
            // sort the adjacent edges ij of i in descending order of x_ij
            auto compare = [&edge_values, &i] (size_t j, size_t k)
            {
                if (i == j)
                    return false;  // i comes last
                if (i == k)
                    return true;  // i comes last
                return edge_values(i, j) > edge_values(i, k);
            };
            std::stable_sort(neighbors.begin(), neighbors.end(), compare);

            for (size_t a = 0; a < n-1; ++a)
            {
                size_t j = neighbors[a];
                if (edge_values(i, j) <= 0.5 + this->min_violation_ / 2)
                    break; // since the neighbors are sorted in descending 
                    // order, there cannot be a violated triangle inequality 
                    // x_ij + x_ik - x_jk <= 1 for any subsequent iteration.
                for (size_t b = a+1; b < n-1; ++b)
                {
                    size_t k = neighbors[b];
                    if (edge_values(i, j) + edge_values(i, k) <= 1 + this->min_violation_)
                        break;  // also because neighbors are sorted
                    double violation = edge_values(i, j) + edge_values(i, k) - edge_values(j, k) - 1;
                    if (violation <= this->min_violation_)
                        continue;  // not enough violation
                    
                    Triangle triangle{{i, j, k}, violation};
                    if (edge2triangle[std::min(i,j)][std::max(i,j)].violation < violation)
                        edge2triangle[std::min(i,j)][std::max(i,j)] = triangle;
                    if (edge2triangle[std::min(i,k)][std::max(i,k)].violation < violation)
                        edge2triangle[std::min(i,k)][std::max(i,k)] = triangle;
                    if (edge2triangle[std::max(j,k)][std::min(j,k)].violation < violation)
                        edge2triangle[std::max(j,k)][std::min(j,k)] = triangle;
                }
            }
        }

        // compute the set of unique violated triangle inequalities
        std::set<Triangle> triangles;
        for (size_t i = 0; i < n; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            if (edge2triangle[i][j].violation > this->min_violation_)
                triangles.insert(edge2triangle[i][j]);
            if (edge2triangle[j][i].violation > this->min_violation_)
                triangles.insert(edge2triangle[j][i]);
        }

        // convert triangle to inequalities
        std::vector<Inequality<int>> inequalities(triangles.size());
        size_t idx = 0;
        for (const Triangle& triangle : triangles)
        {   
            auto t = triangle.triangle;
            inequalities[idx] = {{{t[0], t[1]}, {t[0], t[2]}, {t[1], t[2]}}, {1, 1, -1}, 1, triangle.violation};
            ++idx;
        }
        this->sort_and_reduce_by_violation_depth(inequalities);
        return inequalities;
    }
};


} // namespace CP
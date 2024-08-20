#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cassert>


#pragma once



namespace CP {

/**
 * Class for representing an inequality.
 * The inequality is represented by a vector of edges and a vector of coefficients of equal size.
 * Edges that are not contained in the vector of edges are assumed to have coefficient 0.
 * This allows for a sparse representation of an inequality instead of storing many 0 coefficients.
 * This class also implements some helper functions like evaluating the inequality wrt to
 * some given edge values or computing the dot product of coefficients of two inequalities.
 * The edges and coefficients are always sorted lexicographically in order to make it easier
 * to compare two different inequalities.
 */
template<class C = int>
class Inequality 
{
public:
    typedef std::vector<std::array<size_t, 2>> EDGES;
    typedef std::vector<C> COEFFICIENTS;

    Inequality() {}
    Inequality(const EDGES& edges, const COEFFICIENTS& coefficients, const C& rhs, const double& violation) 
    :
        edges_(edges),
        coefficients_(coefficients),
        rhs_(rhs),
        violation_(violation)
    {
        assert (coefficients_.size() == edges_.size());
        euclidean_ = 0;
        for (C c : coefficients_)
            euclidean_ += c*c;
        euclidean_ = std::sqrt(euclidean_);
        violation_depth_ = violation_ / euclidean_;

        // ensure that source id of edge is smaller than target id
        for (auto& edge : edges_)
        {
            if (edge[0] > edge[1])
            {
                std::swap(edge[0], edge[1]);
            }
            assert (edge[0] < edge[1]);
        }
        // sort edges by ids. For two edges e0 = (s0, t0) and e1 = (s1, t1) we 
        auto compare = [this](size_t i1, size_t i2) {
            if (edges_[i1][0] == edges_[i2][0])
                return edges_[i1][1] < edges_[i2][1];
            else 
                return edges_[i1][0] < edges_[i2][0];
        };
        std::vector<size_t> idx(edges.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), compare);

        EDGES sorted_edges(edges_.size());
        COEFFICIENTS sorted_coefficients(coefficients_.size());
        for (size_t i = 0; i < idx.size(); ++i)
        {
            sorted_edges[i] = edges_[idx[i]];
            sorted_coefficients[i] = coefficients_[idx[i]];
        }
        edges_ = sorted_edges;
        coefficients_ = sorted_coefficients;
    }

    const EDGES& edges() const {return edges_;}
    const COEFFICIENTS& coefficients() const {return coefficients_;}
    const C& rhs() const {return rhs_;}
    const double& violation() const {return violation_;}
    const double& euclidean() const {return euclidean_;}
    const double& violation_depth() const {return violation_depth_;}

    // evaluate the dot product of the coefficients of the inequality and
    // the given edge values
    template<class EDGE_VALUE_MAP>
    double evaluate(const EDGE_VALUE_MAP& edge_values) const
    {
        double value = 0;
        for (size_t i = 0; i < edges_.size(); ++i)
        {
            value += edge_values(edges_[i][0], edges_[i][1]) * coefficients_[i];
        }
        return value;
    }

    // print the inequality to the console
    void print() const
    {
        for (size_t i = 0; i < edges_.size(); ++i)
        {
            std::cout << ((coefficients_[i] > 0) ? " +" : " ") << coefficients_[i] << " x_{" << edges_[i][0] << "," << edges_[i][1] << "}";
        }
        std::cout << " <= " << rhs_ << "\n";
    }

    // compute the dot product of the coefficients of this inequality with the coefficients of another inequality
    C dot_product(const Inequality<C>& other) const
    {
        C dot = 0;
        const EDGES& other_edges = other.edges();
        const COEFFICIENTS& other_coefficients = other.coefficients();
        size_t i = 0;
        size_t j = 0;
        while (i < edges_.size() && j < other_edges.size())
        {
            if (edges_[i][0] == other_edges[j][0])
                if (edges_[i][1] == other_edges[j][1])
                {
                    dot += coefficients_[i] * other_coefficients[j];
                    ++i;
                }
                else if (edges_[i][1] < other_edges[j][1])
                    ++i;
                else
                    ++j;
            else if (edges_[i][0] < other_edges[j][0])
                ++i;
            else
                ++j;
        }
        return dot;
    }


private:

    EDGES edges_;  // edges that have non-zero coefficient in the inequality
    COEFFICIENTS coefficients_;  // coefficients of the edges
    C rhs_;  // right hand side of the inequality

    double violation_;  // amount by which inequality is violated
    double euclidean_; // euclidean distance of inequality
    double violation_depth_; // violation normalized by euclidean distance of inequality

};

} // namespace CP

#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <functional>

#include <inequalities/inequality.hxx>
#include <inequalities/separator.hxx>

#pragma once


namespace CP {


template<class EDGE_VALUE_MAP>
class CompleteSeparator : public AbstractSeparator<int, EDGE_VALUE_MAP> 
{
public:
    CompleteSeparator(size_t n, size_t max_num, size_t k) : 
        AbstractSeparator<int, EDGE_VALUE_MAP>(n, max_num),
        k_(k),
        coefficients_()
    {
        std::ifstream file;

        std::string file_path = "../include/inequalities/complete/cp_facets_" + std::to_string(k_) + ".txt";
        file.open(file_path);

        if (!file.is_open())
            throw std::runtime_error("Unable to open file");

        for (std::string line; std::getline(file, line);)
        {
            coefficients_.push_back({});
            std::istringstream s(line);
            std::copy(
                std::istream_iterator<double>(s), 
                std::istream_iterator<double>(), 
                std::back_inserter(coefficients_.back())  // destination
            );
            assert (coefficients_.back().size() == k * (k-1) / 2 + 1);
        }
    }
    
    std::string name() { return "Complete"+std::to_string(k_); }

    std::vector<Inequality<int>> separate_(EDGE_VALUE_MAP edge_values)
    {
        std::vector<Inequality<int>> inequalities;

        std::vector<size_t> choice(k_);
        iterate_(choice, 0, inequalities, edge_values);

        this->sort_and_reduce_by_euclidean_violation(inequalities);
        this->reduce_by_parallelism(inequalities);

        return inequalities;
    }


private:

    size_t k_;

    std::vector<std::vector<int>> coefficients_;


    void iterate_(std::vector<size_t>& choice, size_t i, std::vector<Inequality<int>>& inequalities, EDGE_VALUE_MAP edge_values)
    {
        size_t start = 0;
        if (i > 0)
            start = choice[i-1] + 1;
        for (choice[i] = start; choice[i] < this->n_; ++choice[i])
        {       
            bool one_exists = false;
            for (size_t j = 0; j < i; ++j)
            { 
                if (edge_values(choice[i], choice[j]) >= 1 - this->min_violation_)
                {
                    one_exists = true;
                    break; 
                }
            }
            if (one_exists)
                continue;
            if (i == 1 && edge_values(choice[0], choice[1]) < this->min_violation_)
                continue; // there must be at least one fractional variable

            if (i == this->k_ - 1)
                separate_choice_(choice, inequalities, edge_values);
            else 
                iterate_(choice, i+1, inequalities, edge_values);
        }
    }

    void separate_choice_(std::vector<size_t>& choice, std::vector<Inequality<int>>& inequalities, EDGE_VALUE_MAP edge_values)
    {
        for (const std::vector<int>& coeff : coefficients_)
        {
            int rhs = coeff.back();
            double lhs = 0;
            for (size_t i = 0; i < k_; ++i)
            for (size_t j = i+1; j < k_; ++j)
            {
                size_t idx = k_ * (k_ - 1) / 2 - (k_ - i)*(k_ - i - 1) / 2 + j - i - 1;
                lhs += edge_values(choice[i], choice[j]) * coeff[idx];
            }
            double violation = lhs - rhs;
            if (violation <= this->min_violation_)
                continue;

            std::vector<std::array<size_t, 2>> edges;
            std::vector<int> edge_coeff;
            for (size_t i = 0; i < k_; ++i)
            for (size_t j = i+1; j < k_; ++j)
            {
                size_t idx = k_ * (k_ - 1) / 2 - (k_ - i)*(k_ - i - 1) / 2 + j - i - 1;
                int c = coeff[idx];
                if (c == 0)
                    continue;
                lhs += edge_values(choice[i], choice[j]) * c;
                edges.push_back({choice[i], choice[j]});
                edge_coeff.push_back(c);
            }
            inequalities.push_back({edges, edge_coeff, rhs, violation});
        }
    }

};


} // namespace CP

#pragma once

#include <vector>
#include <set>

#include <andres/partition.hxx>
#include <kernighan-lin.hxx>


namespace CP 
{

// round fractional solution with kernighan lin based greedy moving algorithm
template<class EDGE_VALUE_MAP, class EDGE_COST_MAP>
typename EDGE_COST_MAP::VALUE_TYPE round_kl(
    const EDGE_VALUE_MAP& edge_values, 
    const EDGE_COST_MAP& edge_costs, 
    size_t steps, 
    std::vector<size_t>& node_labels)
{
    assert (edge_costs.n() == edge_values.n());
    typename EDGE_COST_MAP::VALUE_TYPE best_objective = 0;
    std::vector<size_t> temp_node_labels(edge_costs.n());
    for (size_t step = 0; step < steps; ++step)
    {
        double threshold = 0.5;
        if (steps > 1)
        {
            threshold = 1.0 - (double)step / (steps - 1);
        }
        andres::Partition<size_t> partition(edge_costs.n());
        for (size_t i = 0; i < edge_costs.n(); ++i)
        for (size_t j = i+1; j < edge_costs.n(); ++j)
        {
            if (edge_values(i, j) <= threshold)
                partition.merge(i, j);
        }
        for (size_t i = 0; i < edge_costs.n(); ++i)
            temp_node_labels[i] = partition.find(i);
        typename EDGE_COST_MAP::VALUE_TYPE objective = kernighanLin(edge_costs, temp_node_labels);
        if (objective > best_objective)
        {
            best_objective = objective;
            node_labels = temp_node_labels;
        }
    }
    return best_objective;
}



} // namespace CP
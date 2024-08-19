#pragma once

#include <vector>
#include <set>

#include <andres/partition.hxx>
#include <kernighan-lin.hxx>


namespace CP 
{

template<class EDGE_VALUE_MAP, class EDGE_COST_MAP>
typename EDGE_COST_MAP::VALUE_TYPE round_kl(EDGE_VALUE_MAP edge_values, EDGE_COST_MAP edge_costs, size_t steps)
{
    assert (edge_costs.n() == edge_values.n());
    typename EDGE_COST_MAP::VALUE_TYPE best_objective = 0;
    std::vector<size_t> node_labels(edge_costs.n());
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
            node_labels[i] = partition.find(i);
        typename EDGE_COST_MAP::VALUE_TYPE objective = kernighanLin(edge_costs, node_labels);
        if (objective > best_objective)
            best_objective = objective;
    }
    return best_objective;
}



} // namespace CP
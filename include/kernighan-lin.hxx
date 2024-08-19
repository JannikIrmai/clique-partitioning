#include <vector>
#include <cassert>
#include <limits>
#include <algorithm>

#pragma once

/// copied and modified from https://github.com/bjoern-andres/graph/blob/master/include/andres/graph/multicut/kernighan-lin.hxx

namespace CP {

/**
 * This algorithm is a kernigan lin style greedy moving algorithm for the clique partitioning problem
 * 
*/ 
template<typename EDGE_COST_MAP>
typename EDGE_COST_MAP::VALUE_TYPE kernighanLin(EDGE_COST_MAP edge_costs, std::vector<size_t>& vertex_labels, bool verbose = false)
{
    typedef typename EDGE_COST_MAP::VALUE_TYPE COST_TYPE;

    assert(edge_costs.n() == vertex_labels.size());
    size_t n = edge_costs.n();

    struct Buffers
    {
        Buffers(size_t n) :
            differences(n),
            is_moved(n)
        {}

        std::vector<COST_TYPE> differences;
        std::vector<char> is_moved;
    } buffer(n);

    auto update_bipartition = [&](std::vector<size_t>& A, std::vector<size_t>& B) -> COST_TYPE
    {
        if (A.empty())
            return 0;

        COST_TYPE gain_from_merging = 0;

        // compute differences for set A
        for (size_t i = 0; i < A.size(); ++i)
        {
            COST_TYPE diff = 0;
            buffer.is_moved[A[i]] = 0;
            
            for (auto v : A)
                if (A[i] != v)
                    diff -= edge_costs(A[i], v);

            for (auto v : B)
            {
                diff += edge_costs(A[i], v);
                gain_from_merging += edge_costs(A[i], v);
            }

            buffer.differences[A[i]] = diff;
        }

        // compute differences for set B
        for (size_t i = 0; i < B.size(); ++i)
        {
            COST_TYPE diff = 0;
            buffer.is_moved[B[i]] = 0;

            for (auto v : B)
                if (B[i] != v)
                    diff -= edge_costs(B[i], v);

            for (auto v : A)
                diff += edge_costs(B[i], v);

            buffer.differences[B[i]] = diff;
        }

        struct Move
        {
            int v { -1 };
            COST_TYPE difference { std::numeric_limits<COST_TYPE>::lowest() };
            char new_label;
        };

        COST_TYPE cumulative_diff = 0;
        std::pair<COST_TYPE, size_t> max_move { std::numeric_limits<COST_TYPE>::lowest(), 0 };
        std::vector<Move> moves;
        
        while (true)
        {
            Move m;

            for (auto a : A)
                if (!buffer.is_moved[a] && buffer.differences[a] > m.difference)
                {
                    m.v = a;
                    m.difference = buffer.differences[a];
                    m.new_label = 'B';
                }

            for (auto b : B)
                if (!buffer.is_moved[b] && buffer.differences[b] > m.difference)
                {
                    m.v = b;
                    m.difference = buffer.differences[b];
                    m.new_label = 'A';
                }

            if (m.v == -1)
                break;

            // update differences
            if (m.new_label == 'B')
            {
                for (auto v : A)
                    if (v != m.v && !buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<COST_TYPE>::lowest())
                        buffer.differences[v] += 2.0*edge_costs(v, m.v);

                for (auto v : B)
                    if (!buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<COST_TYPE>::lowest())
                        buffer.differences[v] -= 2.0*edge_costs(v, m.v);
                
            }
            else
            {
                for (auto v : A)
                    if (!buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<COST_TYPE>::lowest())
                        buffer.differences[v] -= 2.0*edge_costs(v, m.v);

                for (auto v : B)
                    if (v != m.v && !buffer.is_moved[v] && buffer.differences[v] > std::numeric_limits<COST_TYPE>::lowest())
                        buffer.differences[v] += 2.0*edge_costs(v, m.v);
            }

            buffer.differences[m.v] = std::numeric_limits<COST_TYPE>::lowest();
            buffer.is_moved[m.v] = 1;
            moves.push_back(m);

            cumulative_diff += m.difference;

            if (cumulative_diff > max_move.first)
                max_move = std::make_pair(cumulative_diff, moves.size());
        }

        if (gain_from_merging > max_move.first && gain_from_merging > 0)
        {
            A.insert(A.end(), B.begin(), B.end());

            B.clear();

            return gain_from_merging;
        }
        if (max_move.first > 0)
        {
            for (size_t i = max_move.second; i < moves.size(); ++i)
                buffer.is_moved[moves[i].v] = 0;

            A.erase(std::partition(A.begin(), A.end(), [&](size_t a) { return !buffer.is_moved[a]; }), A.end());
            B.erase(std::partition(B.begin(), B.end(), [&](size_t b) { return !buffer.is_moved[b]; }), B.end());

            for (size_t i = 0; i < max_move.second; ++i)
                // move vertex to the other set
                if (moves[i].new_label == 'B')
                    B.push_back(moves[i].v);
                else
                    A.push_back(moves[i].v);

            return max_move.first;
        }

        return 0;
    };

    // TODO: Rename decrease to increase
    COST_TYPE starting_cost = 0;

    for(size_t i = 0; i < n; ++i)
    for(size_t j = i+1; j < n; ++j)
    {
        if (vertex_labels[i] == vertex_labels[j])
            starting_cost += edge_costs(i, j);
    }

    size_t numberOfComponents = *std::max_element(vertex_labels.begin(), vertex_labels.end()) + 1;

    // build partitions
    // TODO: If vertex labels are not consecutive this may be problematic
    std::vector<std::vector<size_t>> partitions(numberOfComponents);
    for (size_t i = 0; i < n; ++i)
        partitions[vertex_labels[i]].push_back(i);

    // iteratively update bipartition in order to minimize the total cost of the multicut
    size_t iter = 0;
    while (true)
    {
        ++iter;
        COST_TYPE cost_increase = 0;

        // update pairs of partitions
        for (size_t i = 0; i < partitions.size() - 1; ++i)
            for (auto j = i + 1; j < partitions.size(); ++j)
                if (!partitions[j].empty())
                    cost_increase += update_bipartition(partitions[i], partitions[j]);

        // remove partitions that became empty after the previous step
        auto new_end = std::partition(partitions.begin(), partitions.end(), [](const std::vector<size_t>& s) { return !s.empty(); });
        partitions.resize(new_end - partitions.begin());

        // try to introduce new partitions
        for (size_t i = 0, p_size = partitions.size(); i < p_size; ++i)       
            while (1)
            {
                std::vector<size_t> new_set;
                cost_increase += update_bipartition(partitions[i], new_set);

                if (!new_set.empty())
                    partitions.emplace_back(std::move(new_set));
                else
                    break;
            }

        if (cost_increase == 0)
            break;

        starting_cost += cost_increase;

    }

    for (size_t i = 0; i < partitions.size(); ++i)
        for (size_t j = 0; j < partitions[i].size(); ++j)
            vertex_labels[partitions[i][j]] = i;

    #ifdef DEBUG
        COST_TYPE final_cost = 0;
        for(size_t i = 0; i < n; ++i)
        for(size_t j = i+1; j < n; ++j)
        {
            if (vertex_labels[i] == vertex_labels[j])
                final_cost += edge_costs(i, j);
        }

        if (final_cost != starting_cost)
            throw std::runtime_error("Wrong final cost.");
    #endif

    return starting_cost;

}

} // namespace CP

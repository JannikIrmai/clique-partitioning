#include <vector>
#include <memory>

#include <time.hxx>

#include <inequalities/inequality.hxx>
#include <inequalities/separator.hxx>
#include <inequalities/triangle.hxx>
#include <inequalities/odd-wheel.hxx>
#include <inequalities/odd-bicycle-wheel.hxx>
#include <inequalities/two-partition.hxx>
#include <inequalities/complete.hxx>
#include <inequalities/half-chorded-odd-cycle.hxx>
#include <inequalities/chvatal-gomory.hxx>
#include <inequalities/hypermetric.hxx>


#pragma once


namespace CP {

template<class EPM>
class Callback
{
public:
    Callback(size_t n) : 
        n_(n)
    {}

    std::vector<Inequality<int>> operator()(EPM edge_values, size_t& stage, size_t min_stage = 0)
    {
        std::vector<Inequality<int>> inequalities;
        // iterate over all stages and call the respective separators
        for (stage = 0; stage < stage2idx_.size(); ++stage)
        {
            for (size_t i : stage2idx_[stage])
            {
                std::vector<Inequality<int>> additional_inequalities = separators_[i]->separate(edge_values);;
                inequalities.insert(inequalities.end(), additional_inequalities.begin(), additional_inequalities.end());
            }
            // sort_and_reduce_sorensen(inequalities);  TODO Sort
            if (inequalities.size() > 0 && stage >= min_stage)
                return inequalities;
        }
        return inequalities;
    }

    std::vector<Inequality<int>> operator()(EPM edge_values)
    {
        size_t stage;
        return this->operator()(edge_values, stage);
    }


    void add_separator(std::string name, size_t stage, size_t max_num = 0)
    {
        if (stage > stage2idx_.size())
            throw std::runtime_error("Cannot add separator to stage " + std::to_string(stage) + " since there are no separators in stage " + std::to_string(stage-1));
        if (name == "Triangle")
            separators_.emplace_back(std::make_unique<TriangleSeparator<EPM>>(n_, max_num));
        else if (name == "OddWheel")
            separators_.emplace_back(std::make_unique<OddWheelSeparator<EPM>>(n_, max_num));
        else if (name == "OddBicycleWheel")
            separators_.emplace_back(std::make_unique<OddBicycleWheelSeparator<EPM>>(n_, max_num));
        else if (name == "TwoPartition")
            separators_.emplace_back(std::make_unique<TwoPartitionSeparator<EPM>>(n_, max_num));
        else if (name == "HCO-Cycle")
            separators_.emplace_back(std::make_unique<HalfChordedOddCycleSeparator<EPM>>(n_, max_num));
        else if (name == "Complete4")
            separators_.emplace_back(std::make_unique<CompleteSeparator<EPM>>(n_, max_num, 4));
        else if (name == "Complete5")
            separators_.emplace_back(std::make_unique<CompleteSeparator<EPM>>(n_, max_num, 5));
        else if (name == "ChvatalGomory")
            separators_.emplace_back(std::make_unique<ChvatalGomorySeparator<EPM>>(n_, max_num));
        else if (name == "HyperMetric")
            separators_.emplace_back(std::make_unique<HyperMetricSeparator<EPM>>(n_, max_num));
        else
            throw std::runtime_error("Invalid separator name '" + name + "'.");
        if (stage == stage2idx_.size())
            stage2idx_.push_back({});
        stage2idx_[stage].push_back(separators_.size() - 1);
    }

    template<typename SEPARATOR, typename... Args>
    void add_separator(size_t stage, Args... args)
    {
        if (stage > stage2idx_.size())
            throw std::runtime_error("Cannot add separator to stage " + std::to_string(stage) + " since there are no separators in stage " + std::to_string(stage-1));
        separators_.emplace_back(std::make_unique<SEPARATOR>(n_, args...));
        if (stage == stage2idx_.size())
            stage2idx_.push_back({});
        stage2idx_[stage].push_back(separators_.size() - 1);
    }

    size_t num_separators()
    {
        return separators_.size();
    }

    std::string name(size_t i)
    {
        return separators_[i]->name();
    }

    size_t num_inequalities(size_t i)
    {
        return separators_[i]->total_num();
    }

    size_t num_calls(size_t i)
    {
        return separators_[i]->num_calls();
    }

    double time(size_t i)
    {
        return separators_[i]->total_time();
    }

    size_t num_stages()
    {
        return stage2idx_.size();
    }

private:
    size_t n_;
    std::vector<std::unique_ptr<AbstractSeparator<int, EPM>>> separators_;
    std::vector<std::vector<size_t>> stage2idx_;
};

} // namespace CP
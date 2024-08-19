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

/**
 * This class manages several separation algorithms for the clique partitioning problem.
 * The separation algorithms are arranged in stages: if no algorithm up to stage x finds
 * a violated inequality, only then algorithms from stage x+1 are considered.
 */
template<class EPM>
class SeparatorCallback
{
public:

    // constructor
    SeparatorCallback(size_t n) : 
        n_(n)
    {}

    // this method returns all violated inequalities that were found by the separation algorithms
    std::vector<Inequality<int>> operator()(EPM edge_values, size_t& stage, size_t min_stage = 0)
    {
        std::vector<Inequality<int>> inequalities;
        // iterate over all stages and call the respective separators
        for (stage = 0; stage < stage2idx_.size(); ++stage)
        {
            for (size_t i : stage2idx_[stage])
            {
                // call the i-th separator of the current stage
                std::vector<Inequality<int>> additional_inequalities = separators_[i]->separate(edge_values);;
                // add the inequalities to the list of all inequalities
                inequalities.insert(inequalities.end(), additional_inequalities.begin(), additional_inequalities.end());
            }
            // if in this stage at least one inequality was found and the min_stage is reached, return the inequalities
            if (inequalities.size() > 0 && stage >= min_stage)
                return inequalities;
        }
        return inequalities;
    }

    // overload 
    std::vector<Inequality<int>> operator()(EPM edge_values)
    {
        size_t stage;
        return this->operator()(edge_values, stage);
    }

    // add a separation algorithm based on a string name
    // specify the stage of the separation algorithm and max_num, i.e. the maximum number of inequalities that 
    // should be returned by one call to the separation algorithm
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

    // A more generic method for adding a separation algorithm
    // This method can add any separation algorithm that is implemented as a child class of the AbstractSeparator class
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

    // return the number of separation algorithms
    size_t num_separators()
    {
        return separators_.size();
    }

    // return the name of the i-th separation algorithm
    std::string name(size_t i)
    {
        return separators_[i]->name();
    }

    // return the total number of inequalities that were computed by the i-th algorithm
    size_t num_inequalities(size_t i)
    {
        return separators_[i]->total_num();
    }

    // return the total number of times the i-th algorithm was called
    size_t num_calls(size_t i)
    {
        return separators_[i]->num_calls();
    }

    // return the total runtime of the i-th algorithm
    double time(size_t i)
    {
        return separators_[i]->total_time();
    }

    // return the number of stages
    size_t num_stages()
    {
        return stage2idx_.size();
    }

private:
    size_t n_;  // size of the problem instance
    std::vector<std::unique_ptr<AbstractSeparator<int, EPM>>> separators_;  // separation algorithms
    std::vector<std::vector<size_t>> stage2idx_;  // algorithms ids in each stage
};

} // namespace CP
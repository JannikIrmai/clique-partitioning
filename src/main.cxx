#include <iostream>
#include <clique-partitioning-bnc.hxx>
#include <load.hxx>


typedef CP::EdgePropertyMap<int> EDGE_COST_MAP;
typedef CP::BnC<EDGE_COST_MAP> LP;

int main(int argc, char* argv[])
{
    if (argc != 2)
        throw std::runtime_error("There can only be one argument!");

    std::string instance_path(argv[1]);
    
    EDGE_COST_MAP edge_costs = CP::load_cplib<int>(instance_path);
    LP lp(edge_costs);

    lp.activate_branching = true;
    lp.max_iter_non_basic = 5;
    lp.max_tail_length = 1000;
    lp.tail_threshold = 1;
    lp.separator_callback().add_separator<CP::TriangleSeparator<LP::EPM>>(0);
    // lp.separator_callback().add_separator<CP::OddWheelSeparator<LP::EPM>>(1);
    lp.verbosity = 2;

    lp.add_triangle_inequalities_based_on_negative_cost();
    lp.add_triangle_inequalities_based_on_positive_cost();

    lp.optimize();

    std::cout << "upper bound = " << std::setprecision(10) << lp.bound() << "\n";
    std::cout << "best integer feasible objective = " << lp.best_objective() << "\n";

    std::ofstream results_file("results.txt");
    lp.print_solution(results_file);

    return 0;
}
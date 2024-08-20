#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cassert>

#include <clique-partitioning-bnc.hxx>
#include <property-maps.hxx>


namespace py = pybind11;


typedef CP::EdgePropertyMap<double> ECM;
typedef CP::BnC<ECM> BnC;


std::tuple<double, double, std::vector<double>> bnc(
    size_t n,
    const std::vector<double> edge_costs,
    bool activate_branching,
    size_t max_iter_non_basic,
    size_t max_tail_length,
    double tail_threshold,
    std::vector<std::tuple<std::string, size_t, size_t>> separators,
    bool add_some_triangles,
    size_t verbosity 
){
    assert (edge_costs.size() == n * (n-1) / 2);
    ECM ecm(n, edge_costs);

    BnC bnc(ecm);

    bnc.activate_branching = activate_branching;
    bnc.max_iter_non_basic = max_iter_non_basic;
    bnc.max_tail_length = max_tail_length;
    bnc.tail_threshold = tail_threshold;
    bnc.verbosity = verbosity;

    if (add_some_triangles)
    {
        bnc.add_triangle_inequalities_based_on_negative_cost();
        bnc.add_triangle_inequalities_based_on_positive_cost();
    }

    for (auto separator : separators)
        bnc.separator_callback().add_separator(std::get<0>(separator), std::get<1>(separator), std::get<2>(separator));

    bnc.optimize();

    std::vector<double> edge_values;
    for (size_t i = 0; i < n; ++i)
    for (size_t j = i+1; j < n; ++j)
        edge_values.push_back(bnc.get_solution(i, j, activate_branching));
    return {bnc.best_objective(), bnc.bound(), edge_values};
}


// wrap as Python module
PYBIND11_MODULE(clique_partitioning, m)
{
    m.def("bnc", &bnc, "Branch and cut algorithm for clique partitioning",
        py::arg("n"),
        py::arg("costs"), 
        py::arg("activate_branching")=true,
        py::arg("max_iter_non_basic")=5,
        py::arg("max_tail_length")=4,
        py::arg("tail_threshold")=0.999,
        py::arg("separators")=std::vector<std::tuple<std::string, size_t, size_t>>({{"Triangle", 0, 0}}),
        py::arg("add_some_triangles")=true,
        py::arg("verbosity")=2); 
}
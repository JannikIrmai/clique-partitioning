#include <inequalities/half-chorded-odd-cycle.hxx>
#include <property-maps.hxx>


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


void test_half_chorded_odd_cycle()
{
    typedef CP::EdgePropertyMap<double> EPM;

    std::vector<double> edge_values = {0.5, 1.0, 1.0, 0.5, 0.5, 1.0, 1.0, 0.5, 1.0, 0.5};
    EPM edge_value_map(5, edge_values);

    CP::HalfChordedOddCycleSeparator<EPM> separator(5);
    std::vector<CP::Cycle> cycles = separator.separate_cycles(edge_value_map);

    test (cycles.size() == 5);
    test (cycles[0].violation == 0.5);
    test (cycles[0].cycle.size() == 5);
    test (cycles[0].cycle[0] == 1);
    test (cycles[0].cycle[1] == 0);
    test (cycles[0].cycle[2] == 4);
    test (cycles[0].cycle[3] == 3);
    test (cycles[0].cycle[4] == 2);

    std::vector<CP::Inequality<int>> inequalities = separator.separate(edge_value_map);
    test (inequalities.size() == 1);
    test (inequalities[0].violation() == 0.5);
    test (inequalities[0].euclidean_violation() == 0.5 / std::sqrt(10));
    test (inequalities[0].evaluate(edge_value_map) == 2.5);

}


int main()
{

    test_half_chorded_odd_cycle();

    return 0;
}
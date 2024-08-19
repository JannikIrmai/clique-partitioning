#include <inequalities/complete.hxx>
#include <property-maps.hxx>
#include <clique-partitioning-bnc.hxx>


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


typedef CP::EdgePropertyMap<std::vector<double>> EPM;


void test_constructor()
{   
    size_t n = 10;
    CP::CompleteSeparator<EPM> sep(n, 0, 5);

    test(sep.name() == "Complete5");

    std::vector<double> edge_values(n*(n-1)/2, 0);
    EPM epm(n, edge_values);
    auto ineq = sep.separate_(epm);
    assert (ineq.size() == 0);
}


void test_bow_tie()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;

    std::vector<double> edge_values = {0.5, 0.5, 0.5, 0.5, 1, 0, 0, 0, 0, 1};
    EPM edge_value_map(5, edge_values);

    CP::CompleteSeparator<EPM> separator(5, 0, 5);
    std::vector<CP::Inequality<int>> inequalities = separator.separate_(edge_value_map);

    // test (inequalities.size() == 1);
}

void test_soup_fractional()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;

    size_t n = 11;
    std::vector<double> edge_values(n * (n - 1) / 2);
    EPM edge_value_map(n, edge_values);
    edge_value_map(0, 1) = 0.5;
    edge_value_map(1, 2) = 0.5;
    edge_value_map(1, 3) = 0.5;
    edge_value_map(2, 3) = 0.5;
    edge_value_map(3, 4) = 0.5;
    edge_value_map(4, 5) = 0.5;
    edge_value_map(4, 6) = 0.5;
    edge_value_map(5, 6) = 0.5;
    edge_value_map(5, 7) = 0.5;
    edge_value_map(5, 9) = 0.5;
    edge_value_map(6, 7) = 0.5;
    edge_value_map(6, 8) = 0.5;
    edge_value_map(7, 8) = 0.5;
    edge_value_map(7, 9) = 0.5;
    edge_value_map(9, 10) = 0.5;

    CP::CompleteSeparator<EPM> separator(n, 0, 5);
    std::vector<CP::Inequality<int>> inequalities = separator.separate_(edge_value_map);

    for (const auto& ineq : inequalities)
        ineq.print();
    
}


int main()
{  
    test_constructor();
    test_bow_tie();
    test_soup_fractional();

    return 0;
}
#include <callback.hxx>
#include <property-maps.hxx>
#include <cassert>


void test_callback()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;
    std::vector<double> edge_values = {0.5, 0.5, 0.5, 0, 0, 0};
    EPM edge_value_map(4, edge_values);

    CP::Callback<EPM> callback(4);
    callback.add_separator("Triangle", 0);
    callback.add_separator("OddWheel", 1);

    auto inequalities = callback(edge_value_map);

    assert (inequalities.size() == 1);
}


int main()
{

    test_callback();
    return 0;
}
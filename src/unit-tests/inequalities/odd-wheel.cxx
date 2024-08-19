#include <inequalities/odd-wheel.hxx>
#include <property-maps.hxx>



void test_3_wheel()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;

    std::vector<double> edge_values = {0.5, 0.5, 0.5, 0, 0, 0};
    EPM edge_value_map(4, edge_values);

    CP::OddWheelSeparator<EPM> separator(4);
    std::vector<CP::Wheel> wheels = separator.separate_wheel(edge_value_map);

    assert (wheels.size() == 1);
    assert (wheels[0].center == 0);
    assert (wheels[0].violation == 0.5);
    assert (wheels[0].wheel[0] == 1);
    assert (wheels[0].wheel[1] == 2);
    assert (wheels[0].wheel[2] == 3);
}


int main()
{

    test_3_wheel();

    return 0;
}
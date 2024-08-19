#include <property-maps.hxx>
#include <vector>
#include <cassert>


void test_edge_property_map()
{
    typedef std::vector<double> VECTOR;
    VECTOR data{1.0, 1.2, 2.4, -4.2, 3.2, 3.3};

    CP::EdgePropertyMap<double> epm(4, data);

    assert (epm(0, 1) == data[0]);
    assert (epm(0, 2) == data[1]);
    assert (epm(0, 3) == data[2]);
    assert (epm(1, 2) == data[3]);
    assert (epm(1, 3) == data[4]);
    assert (epm(2, 3) == data[5]);

    epm(1, 3) = 8;
    assert(data[4] == 3.2);
}

int main()
{
    test_edge_property_map();

    return 0;
}

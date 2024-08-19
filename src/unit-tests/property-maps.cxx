#include <property-maps.hxx>
#include <vector>
#include <cassert>


void test_edge_property_map_ref()
{
    typedef std::vector<double> VECTOR;
    VECTOR vector{1.0, 1.2, 2.4, -4.2, 3.2, 3.3};

    CP::EdgePropertyMap<VECTOR> epm(4, vector);

    assert (epm(0, 1) == vector[0]);
    assert (epm(0, 2) == vector[1]);
    assert (epm(0, 3) == vector[2]);
    assert (epm(1, 2) == vector[3]);
    assert (epm(1, 3) == vector[4]);
    assert (epm(2, 3) == vector[5]);

    epm(1, 3) = 8;
    assert(vector[4] == 8);
}

int main()
{
    test_edge_property_map_ref();

    return 0;
}

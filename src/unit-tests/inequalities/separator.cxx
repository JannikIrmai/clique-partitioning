#include <vector>
#include <memory>
#include <cassert>

#include <property-maps.hxx>
#include <inequalities/separator.hxx>
#include <inequalities/triangle.hxx>
#include <inequalities/odd-wheel.hxx>


void test_vector_of_separators()
{
    typedef CP::EdgePropertyMap<std::vector<double>> EPM;
    typedef CP::AbstractSeparator<int, EPM> AbstractSeparator;
    typedef CP::TriangleSeparator<EPM> TriangleSeparator;
    typedef CP::OddWheelSeparator<EPM> OddWheelSeparator;

    size_t n = 3;
    std::vector<std::unique_ptr<AbstractSeparator>> vector;

    vector.emplace_back(std::make_unique<TriangleSeparator>(n));
    vector.emplace_back(std::make_unique<OddWheelSeparator>(n));

    assert (vector[0]->name() == "Triangle");
    assert (vector[1]->name() == "OddWheel");   
}


int main()
{
    test_vector_of_separators();

    return 0;
}
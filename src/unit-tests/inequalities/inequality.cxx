#include <cassert>
#include <inequalities/inequality.hxx>



void test_dot_product()
{
    typedef int C;
    typedef std::vector<std::array<size_t, 2>> EDGES;
    typedef std::vector<C> COEFFICIENTS;

    EDGES edges0 = {{0, 1}, {1, 2}, {2, 3}};
    COEFFICIENTS coefficients0 = {1, 2, 3};

    EDGES edges1 = {{0, 3}, {3, 2}, {0, 1}, {0, 2}, {3, 1}};
    COEFFICIENTS coefficients1 = {-1, -2, -3, 1, 2};

    CP::Inequality<C> inequality0(edges0, coefficients0, 0, 0);
    CP::Inequality<C> inequality1(edges1, coefficients1, 0, 0);

    assert (inequality0.dot_product(inequality1) == 1 * (-3) + 3 * (-2));
}



int main()
{
    test_dot_product();

    return 0;
}
#include <load.hxx>
#include <cassert>


void test_load()
{
    auto data = CP::load_cplib<double>("../data/example-from-paper.txt");

    assert (data.n() == 6);
    assert (data(0, 1) == 1);
    assert (data(0, 2) == 3);
    assert (data(4, 5) == -2);
}


int main()
{
    test_load();
    return 0;
}
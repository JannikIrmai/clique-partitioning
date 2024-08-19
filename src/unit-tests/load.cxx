#include <load.hxx>
#include <cassert>


void test_load()
{
    auto data = CP::load_cplib<double>("../data/example-from-paper.txt");

    assert (data.first == 6);
    assert (data.second.size() == 15);
}


int main()
{
    test_load();
    return 0;
}
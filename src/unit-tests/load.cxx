#include <load.hxx>
#include <cassert>


void test_load()
{
    auto data = CP::load_cplib<double>("../Modularity/karate.txt");

    assert (data.first == 34);
    assert (data.second.size() == 561);
}


int main()
{
    test_load();
    return 0;
}
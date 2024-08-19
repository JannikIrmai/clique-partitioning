#include <inequalities/two-partition.hxx>
#include <property-maps.hxx>


void test(const bool& b)
{
    if (!b)
        throw std::logic_error("Test Failed!");
}


void test_two_partition()
{
    typedef CP::EdgePropertyMap<double> EPM;

    std::vector<double> edge_values = {0.5, 0.5, 0.5, 0, 0, 0};
    EPM edge_value_map(4, edge_values);

    CP::TwoPartitionSeparator<EPM> separator(4);
    std::vector<CP::TwoPartition> two_partitions = separator.separate_two_partition(edge_value_map);

    auto ineq = separator.separate(edge_value_map);

    test (two_partitions.size() == 3);
    test (two_partitions[0].s.size() == 1);
    test (two_partitions[0].t.size() == 3);
    test (two_partitions[0].violation == 0.5);
    test (two_partitions[0].s[0] == 0);
    test (two_partitions[0].t[0] == 1);
    test (two_partitions[0].t[1] == 2);
    test (two_partitions[0].t[2] == 3);
}


int main()
{
    test_two_partition();

    return 0;
}
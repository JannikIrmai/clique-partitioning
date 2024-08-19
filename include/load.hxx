#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <property-maps.hxx>


namespace CP {

// load a CP-Lib instance in the format of https://github.com/MMSorensen/CP-Lib
template<class T>
EdgePropertyMap<T> load_cplib(const std::string file_name)
{
    // open text file
    std::ifstream file(file_name);
    std::string line;
    
    if (!file.is_open())
        throw std::runtime_error("Unable to open file");
    
    // read number of nodes
    std::getline(file, line);
    size_t n = std::stoul(line);

    // read edge costs to vector
    std::vector<T> edge_costs;

    for (size_t i = 0; i < n-1; ++i)
    {
        std::getline(file, line);
        std::istringstream s(line);
        std::copy(
            std::istream_iterator<double>(s), 
            std::istream_iterator<double>(), 
            std::back_inserter(edge_costs)
        );
    }

    // return number of nodes and edges costs
    return {n, edge_costs};
}

} // namespace load


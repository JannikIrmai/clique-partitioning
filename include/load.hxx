#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>


namespace CP {

// load a CP-Lib instance in the format of https://github.com/MMSorensen/CP-Lib

template<class T>
std::pair<size_t, std::vector<T>> load_cplib(const std::string file_name)
{

    std::ifstream file;
    std::string line;

    file.open(file_name);

    if (!file.is_open())
        throw std::runtime_error("Unable to open file");
    
    // read number of nodes
    std::getline(file, line);
    size_t n = std::stoul(line);

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

    return {n, edge_costs};
}

} // namespace load


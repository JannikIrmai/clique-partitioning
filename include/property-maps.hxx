#pragma once

#include <cstddef>
#include <cassert>
#include <utility>
#include <vector>

namespace CP 
{

/**
 * This class implements a data structure for storing properties of edges of a complete graph with n node.
 * The properties are stored in a continuous vector of length n*(n-1)/2.
 * The operator() allows read and write access of these elements.
 */
template<typename T>
class EdgePropertyMap
{
public:

    typedef T VALUE_TYPE;

    EdgePropertyMap(size_t n, const T& default_value) :
        n_(n),
        data_(n*(n-1)/2, default_value)
    {}

    EdgePropertyMap(size_t n) :
        EdgePropertyMap(n, T())
    {}

    template<class DATA>
    EdgePropertyMap(size_t n, const DATA& data) :
        n_(n),
        data_(data)
    {
        assert (data.size() == n*(n-1)/2);
    }

    VALUE_TYPE operator()(size_t i, size_t j) const
    {
        size_t idx = edge2idx_(i, j);
        return data_[idx];
    }

    VALUE_TYPE& operator()(size_t i, size_t j)
    {
        size_t idx = edge2idx_(i, j);
        return data_[idx];
    }

    size_t n() const
    {
        return n_;
    }

private:
    const size_t n_;
    std::vector<VALUE_TYPE> data_;

    size_t edge2idx_(size_t i, size_t j) const
    {
        assert (i != j);
        if (i > j)
            std::swap(i, j);
        assert (j < n_);
        return data_.size() - (n_-i)*(n_-i-1) / 2 + j - i - 1;
    }

};

} // namespace CP
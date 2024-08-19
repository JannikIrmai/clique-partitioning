#pragma once

#include <cstddef>
#include <cassert>
#include <utility>
#include <vector>

namespace CP 
{

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

    template<class OTHER>
    EdgePropertyMap(const OTHER& other) :
        n_(other.n()),
        data_(other.data())
    {}

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

    const std::vector<VALUE_TYPE>& data() const
    {
        return data;
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
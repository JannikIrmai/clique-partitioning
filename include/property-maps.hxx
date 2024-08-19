#include <cstddef>
#include <cassert>
#include <utility>

#pragma once


namespace CP 
{

template<class RANDOM_ACCESS_CONTAINER>
class EdgePropertyMap
{
public:
    typedef typename RANDOM_ACCESS_CONTAINER::value_type VALUE_TYPE; 

    EdgePropertyMap( size_t n, RANDOM_ACCESS_CONTAINER& data) :
        n_(n),
        data_(data)
    {
        assert (data.size() == (n * (n-1)) / 2);
    }

    VALUE_TYPE& operator()(size_t i, size_t j) const
    {
        assert (i != j);
        if (i > j)
            std::swap(i, j);
        assert (j < n_);

        size_t idx = data_.size() - (n_-i)*(n_-i-1) / 2 + j - i - 1;
        return data_[idx];
    }

    size_t n()
    {
        return n_;
    }

private:
    const size_t n_;
    RANDOM_ACCESS_CONTAINER& data_;

};



} // namespace CP
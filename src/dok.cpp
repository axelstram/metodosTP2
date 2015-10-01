#include "dok.h"

DictionaryOfKeys::DictionaryOfKeys(size_t n, size_t m) : rows_(n), cols_(m)
{
}


double& DictionaryOfKeys::operator()(size_t i, size_t j)
{
	return data_[make_pair(i,j)]; // http://www.cplusplus.com/reference/map/map/operator%5B%5D/
}



double DictionaryOfKeys::operator()(size_t i, size_t j) const
{
    return data_.count(make_pair(i,j)) == 1 ? data_.find(make_pair(i,j))->second : 0.0;
} 


size_t DictionaryOfKeys::rows() const
{
	return rows_;
}


size_t DictionaryOfKeys::cols() const
{
	return cols_;
}

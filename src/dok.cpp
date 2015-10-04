#include "dok.h"
#include <cassert>

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


DictionaryOfKeys* DictionaryOfKeys::operator*(double scalar)
{
	DictionaryOfKeys& thisMat = *this;
	DictionaryOfKeys* res = new DictionaryOfKeys(thisMat.cols(), thisMat.rows());
	
	if(scalar != 0){
  		for (map<pair<int,int>,double>::iterator it = data_.begin(); it != data_.end(); ++it)
			res->operator()(it->first.first,it->first.second) = it->second*scalar;
	}

	return res;
}

vector<double> DictionaryOfKeys::operator*(const vector<double>& x)
{
	vector<double> res;
	DictionaryOfKeys& thisMatrix = *this;

	assert(x.size() == thisMatrix.cols() && "Error al multiplicar matriz y vector de diferentes dimensiones");

	for (int i = 0; i < thisMatrix.rows(); i++) {
		double mult = 0;
		for(int j = 0; j < thisMatrix.cols(); j++) {
			mult += thisMatrix(i,j) * x[j];
		}
		res.push_back(mult);
	}

	return res;
}
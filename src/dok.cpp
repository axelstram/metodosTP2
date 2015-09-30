#include "dok.h"

DictionaryOfKeys::DictionaryOfKeys(size_t n, size_t m) : rows_(n), cols_(m)
{
}


DictionaryOfKeys::~DictionaryOfKeys()
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


/*
DictionaryOfKeys& DictionaryOfKeys::operator+(const DictionaryOfKeys& anotherDictionaryOfKeys)
{

}



DictionaryOfKeys& DictionaryOfKeys::operator-(const DictionaryOfKeys& anotherDictionaryOfKeys)
{

}



DictionaryOfKeys& DictionaryOfKeys::operator*(const DictionaryOfKeys& anotherDictionaryOfKeys)
{
	//assert(this->cols_ == anotherDictionaryOfKeys.rows_ && "Cols and Rows don't DictionaryOfKeysch");
}

*/

size_t DictionaryOfKeys::rows() const
{
	return rows_;
}



size_t DictionaryOfKeys::cols() const
{
	return cols_;
}


/*
DictionaryOfKeys DictionaryOfKeys::clone() const
{
	DictionaryOfKeys res(rows_, cols_);
	const DictionaryOfKeys& thisDictionaryOfKeys = *this;

	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			res(i, j) = thisDictionaryOfKeys(i, j);
		}
	}

	return res;
}


*/

void DictionaryOfKeys::ShowOctave()
{
	DictionaryOfKeys& thisDictionaryOfKeys = *this;
	cout << "[";
	for (int j = 0; j < rows_; j++) {
		for (int k = 0; k < cols_; k++) {
			if (thisDictionaryOfKeys(j,k) == 1 || thisDictionaryOfKeys(j,k) == 0)
				cout << thisDictionaryOfKeys(j,k) << ".000000000 ";
			else
				cout << thisDictionaryOfKeys(j,k) << " ";
		}
		cout << "; ";
	}
	cout << "]";
}



void DictionaryOfKeys::Show()
{
	DictionaryOfKeys& thisDictionaryOfKeys = *this;
	for (int j = 0; j < rows_; j++) {
		for (int k = 0; k < cols_; k++) {
			if (thisDictionaryOfKeys(j,k) == 1 || thisDictionaryOfKeys(j,k) == 0)
				cout << thisDictionaryOfKeys(j,k) << ".000000000 ";
			else
				cout << thisDictionaryOfKeys(j,k) << " ";
		}
		cout << endl;
	}
}
/*
void DictionaryOfKeys::GetKeys(){
std::vector<char> v;
	for(pair<int,int>,double>::iterator it = data_.begin(); it != data_.end(); ++it) {
		v.push_back(it->first);
		cout << it->first << "\n";
	}
}
*/

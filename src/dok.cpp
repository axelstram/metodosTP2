#include "DictionaryOfKeys.h"


DictionaryOfKeys::DictionaryOfKeys(std::vector<pair< ,std::vector<int> v;> > v;)
{



}



double& DictionaryOfKeys::operator()(size_t i, size_t j)
{
    return 0;
}



double DictionaryOfKeys::operator()(size_t i, size_t j) const
{
    return 0;
} 



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



size_t DictionaryOfKeys::rows() const
{
	return rows_;
}



size_t DictionaryOfKeys::cols() const
{
	return cols_;
}



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
		cout << ";";
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
// Redefinir operador <<

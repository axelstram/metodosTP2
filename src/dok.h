#include <iostream>
#include <map> 
#include <utility> // pair
#include "matrix.h" 


using namespace std;
class DictionaryOfKeys : public Matrix
{
	public:
		DictionaryOfKeys(size_t n, size_t m);

		double& operator()(size_t i, size_t j);
		double operator()(size_t i, size_t j) const;
		size_t rows() const;
		size_t cols() const;
		void Show();
		void ShowOctave();

	private:
		size_t rows_;
		size_t cols_;
		map<pair<int,int>,double> data_;
};


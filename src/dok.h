#include <iostream>
#include <map>
#include <utility> // pair
// #include <vector>
// #include <cassert>
// #include <fstream>


using namespace std;
class DictionaryOfKeys
{
	public:
	DictionaryOfKeys(size_t n, size_t m);
	~DictionaryOfKeys();

	double& operator()(size_t i, size_t j);
	double operator()(size_t i, size_t j) const;
//	Mat& operator+(const Mat& anotherMat);
//	Mat& operator-(const Mat& anotherMat);
//	Mat& operator*(const Mat& anotherMat);
	size_t rows() const;
	size_t cols() const;
	//Mat clone() const;

	void Show();
	void ShowOctave();

	private:

	size_t rows_;
	size_t cols_;
	map<pair<int,int>,double> data_;
};


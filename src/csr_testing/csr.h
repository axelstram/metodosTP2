#include <iostream>
#include <vector>
// #include <cassert>
// #include <fstream>

using namespace std;


class CompressedSparseRow
{
	public:
	CompressedSparseRow(size_t n, size_t m);
	~CompressedSparseRow();

	double operator()(size_t i, size_t j);
	void operator()(size_t i, size_t j, double value);
//	Mat& operator+(const Mat& anotherMat);
//	Mat& operator-(const Mat& anotherMat);
//	Mat& operator*(const Mat& anotherMat);
	size_t rows() const;
	size_t cols() const;
	//Mat clone() const;

	void Show();
	void show_vectors();
	// void ShowOctave();

	private:
	size_t rows_; 
	size_t cols_;

	vector<double> value_; //non-zero values
	vector<double> row_ptr_; //locations in the val vector that start a row
	vector<int> col_ind_; //colum indexes
};


#include <iostream>
#include <vector>
#include "matrix.h"	
using namespace std;


class CompressedSparseRow : public Matrix
{
	public:
		CompressedSparseRow(size_t n, size_t m);
		//~CompressedSparseRow();

		double operator()(size_t i, size_t j) const; //leer
		double & operator()(size_t i, size_t j); //agregar elemento 
		size_t rows() const;
		size_t cols() const;
		void Show();
		void ShowOctave();

		void show_vectors();
		//Mat& operator+(const Mat& anotherMat);
		//Mat& operator-(const Mat& anotherMat);
		//Mat& operator*(const Mat& anotherMat);
		//Mat clone() const;
	private:
	size_t rows_;
	size_t cols_;

	vector<double> value_; //non-zero values
	vector<int> row_ptr_; //locations in the val vector that start a row
	vector<int> col_ind_; //colum indexes
};


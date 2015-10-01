/***** A 2-D Matrix *****/
#ifndef MAT_H_
#define MAT_H_

#include "matrix.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

class Mat : public Matrix{
	public:
	    Mat(size_t rows_, size_t cols_);
	    double& operator()(size_t i, size_t j);
	    double operator()(size_t i, size_t j) const;
	    Mat& operator+(const Mat& anotherMat);
	    Mat& operator-(const Mat& anotherMat);
	    Mat& operator*(const Mat& anotherMat);
	    Mat& operator*(double scalar);
	    vector<double> operator*(const vector<double>& x);
	    size_t rows() const;
	    size_t cols() const;
	    Mat clone() const;
		void Show();
		void ShowOctave();

	private:
	    size_t rows_;
	    size_t cols_;
	    vector<double> data_;
	};
 




#endif

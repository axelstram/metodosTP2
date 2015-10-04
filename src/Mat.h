/***** A 2-D Matrix *****/
#ifndef MAT_H_
#define MAT_H_


#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include "matrix.h"

using namespace std;


class Mat : public Matrix{
	public:
	    Mat(size_t rows_, size_t cols_);
	    Mat(const Mat& anotherMat);
	    double& operator()(size_t i, size_t j);
	    double operator()(size_t i, size_t j) const;
	    vector<double> operator*(const vector<double>& x);
	    Mat* operator*(double scalar);
	    size_t rows() const;
	    size_t cols() const;
		void Show();
		void ShowOctave();
	    //Mat operator+(const Mat& anotherMat);
	    //Mat operator*(const Mat& anotherMat);
	    //Mat& operator=(const Mat& anotherMat);
	    //Mat clone() const;
		//void checkBandMat();

	private:

	    size_t rows_;
	    size_t cols_;
	    vector<double> data_;
	};
 




#endif

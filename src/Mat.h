/***** A 2-D Matrix *****/
#ifndef MATRIX_H_
#define MATRIX_H_


#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>


using namespace std;


class Mat {
	public:
	    Mat(size_t rows_, size_t cols_);
	    double& operator()(size_t i, size_t j);
	    double operator()(size_t i, size_t j) const;
	    Mat& operator+(const Mat& anotherMat);
	    Mat& operator-(const Mat& anotherMat);
	    Mat& operator*(const Mat& anotherMat);
	    size_t rows() const;
	    size_t cols() const;
	    Mat clone() const;

		void Show();
		void ShowOctave();
		void checkBandMat();

	private:
		bool is_band;
		int band_ceros_count;

	    size_t rows_;
	    size_t cols_;
	    vector<double> data_;
	};
 




#endif

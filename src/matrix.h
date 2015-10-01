/***** A 2-D Matrix *****/
#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>

using namespace std;

class Matrix {
	protected:
    	size_t rows_;
		size_t cols_;
	
	public:
	    virtual double& operator()(size_t i, size_t j) =0; //metodo puramente virtual
	    virtual double operator()(size_t i, size_t j) const =0;
	    vector<double> operator*(const vector<double>& x);
	    virtual size_t rows() const =0; 
	    virtual size_t cols() const =0;
		void Show();
		void ShowOctave();
	};

#endif

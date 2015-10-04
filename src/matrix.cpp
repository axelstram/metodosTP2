#include "matrix.h"
#include <cassert>
#include <iostream>

/*
vector<double> Matrix::operator*(const vector<double>& x)
{
	vector<double> res;
	Matrix& thisMatrix = *this;

	assert(x.size() == thisMatrix.cols() && "Error al multiplicar matriz y vector de diferentes dimensiones");

	for (int i = 0; i < thisMatrix.rows(); i++) {
		double mult = 0;
		for(int j = 0; j < thisMatrix.cols(); j++) {
			mult += thisMatrix(i,j) * x[j];
		}
		res.push_back(mult);
	}

	return res;
}
*/

void Matrix::ShowOctave()
{
	Matrix& thisMatrix = *this;
	cout << "[";
	for (int j = 0; j < thisMatrix.rows(); j++) {
		for (int k = 0; k < thisMatrix.cols(); k++) {
			if (thisMatrix(j,k) == 1 || thisMatrix(j,k) == 0)
				cout << thisMatrix(j,k) << ".000000000 ";
			else
				cout << thisMatrix(j,k) << " ";
		}
		cout << "; ";
	}
	cout << "]";
}



void Matrix::Show()
{
	Matrix& thisMatrix = *this;
	for (int j = 0; j < thisMatrix.rows(); j++) {
		for (int k = 0; k < thisMatrix.cols(); k++) {
			if (thisMatrix(j,k) == 1 || thisMatrix(j,k) == 0)
				cout << thisMatrix(j,k) << ".000000000 ";
			else
				cout << thisMatrix(j,k) << " ";
		}
		cout << endl;
	}
}
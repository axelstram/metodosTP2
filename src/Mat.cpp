#include "Mat.h"


Mat::Mat(size_t n, size_t m) : rows_(n), cols_(m), data_((n) * (m))
{
}



double& Mat::operator()(size_t i, size_t j)
{
    return data_[i * cols_ + j];
}



double Mat::operator()(size_t i, size_t j) const
{
    return data_[i * cols_ + j];
} 



Mat Mat::operator+(const Mat& anotherMat)
{
	Mat& thisMat = *this;

	Mat res(rows_, cols_);

	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			res(i, j) = thisMat(i, j) + anotherMat(i, j);
		}
	}

	return res;
}



Mat::Mat(const Mat& anotherMat)
{
	this->rows_ = anotherMat.rows_;
	this->cols_ = anotherMat.cols_;
	this->data_ = anotherMat.data_;
}



Mat& Mat::operator=(const Mat& anotherMat)
{
	this->rows_ = anotherMat.rows_;
	this->cols_ = anotherMat.cols_;
	this->data_ = anotherMat.data_;

	return *this;
}



Mat Mat::operator*(const Mat& anotherMat)
{
	//assert(this->cols_ == anotherMat.rows_ && "Cols and Rows don't match");
}



Mat Mat::operator*(double scalar)
{
	Mat& thisMat = *this;
	Mat res(cols_, rows_);

	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			res(i, j) = thisMat(i, j) * scalar;
		}
	}

	return res;
}


vector<double> Mat::operator*(const vector<double>& x)
{
	vector<double> res(x.size());
	Mat& thisMat = *this;
	
	assert(x.size() == cols_&& "Error al multiplicar matriz y vector de diferentes dimensiones");

	for (int i = 0; i < rows_; i++) {
		double mult = 0;
		for(int j = 0; j < cols_; j++) {
			mult += thisMat(i,j) * x[j];
		}
		res[i] = mult;
	}

	return res;
}


size_t Mat::rows() const
{
	return rows_;
}



size_t Mat::cols() const
{
	return cols_;
}



Mat Mat::clone() const
{
	Mat res(rows_, cols_);
	const Mat& thisMat = *this;

	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			res(i, j) = thisMat(i, j);
		}
	}

	return res;
}



void Mat::ShowOctave()
{
	Mat& thisMat = *this;
	cout << "[";
	for (int j = 0; j < rows_; j++) {
		for (int k = 0; k < cols_; k++) {
			if (thisMat(j,k) == 1 || thisMat(j,k) == 0)
				cout << thisMat(j,k) << ".000000000 ";
			else
				cout << thisMat(j,k) << " ";
		}
		cout << ";";
	}
	cout << "]";
}



void Mat::Show()
{


	Mat& thisMat = *this;
	for (int j = 0; j < rows_; j++) {
		for (int k = 0; k < cols_; k++) {
			if (thisMat(j,k) == 1 || thisMat(j,k) == 0)
				cout << thisMat(j,k) << ".0 ";
				//cout << thisMat(j,k) << ".000000000 ";
			else
				cout << thisMat(j,k) << " ";
		}
		cout << endl;
	}

}
// Redefinir operador <<

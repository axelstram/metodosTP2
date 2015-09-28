#include "csr.h"


CompressedSparseRow::CompressedSparseRow(size_t n, size_t m) : rows_(n), cols_(m) {}


CompressedSparseRow::~CompressedSparseRow() {}


//usado para leer un elemento
double CompressedSparseRow::operator()(size_t i, size_t j)
{
	int idx_init_row = row_ptr_[i];
	//busco la posicion de la columna j en el area de la fila i
	for(int k = idx_init_row; k < row_ptr_[i+1]-1; k++){
		if(j == col_ind_[k]) return value_[k];
	}
	return 0;
}


//para agregar un elemento
void CompressedSparseRow::operator()(size_t i, size_t j, double value)
{
	//caso en que hay una o varias filas de ceros
	if(col_ind_.size()-1 < i){
		for(int k = 0; k < i-col_ind_.size(); k++){
			col_ind_.push_back(-1);
			value_.push_back(-1.0);
			row_ptr_.push_back(value_.size());
		}
	}    
	//add add value an its column at the end (for initialization only)
	col_ind_.push_back(j);
	value_.push_back(value);
	row_ptr_.push_back(value_.size()-1);

} 


/*
CompressedSparseRow& CompressedSparseRow::operator+(const CompressedSparseRow& anotherCompressedSparseRow)
{

}



CompressedSparseRow& CompressedSparseRow::operator-(const CompressedSparseRow& anotherCompressedSparseRow)
{

}



CompressedSparseRow& CompressedSparseRow::operator*(const CompressedSparseRow& anotherCompressedSparseRow)
{
	//assert(this->cols_ == anotherCompressedSparseRow.rows_ && "Cols and Rows don't CompressedSparseRowch");
}

*/

size_t CompressedSparseRow::rows() const
{
	return rows_;
}



size_t CompressedSparseRow::cols() const
{
	return cols_;
}


/*
CompressedSparseRow CompressedSparseRow::clone() const
{
	CompressedSparseRow res(rows_, cols_);
	const CompressedSparseRow& thisCompressedSparseRow = *this;

	for (int i = 0; i < rows_; i++) {
		for (int j = 0; j < cols_; j++) {
			res(i, j) = thisCompressedSparseRow(i, j);
		}
	}

	return res;
}


*/

// void CompressedSparseRow::ShowOctave()
// {
// 	CompressedSparseRow& thisCompressedSparseRow = *this;
// 	cout << "[";
// 	for (int j = 0; j < rows_; j++) {
// 		for (int k = 0; k < cols_; k++) {
// 			if (thisCompressedSparseRow(j,k) == 1 || thisCompressedSparseRow(j,k) == 0)
// 				cout << thisCompressedSparseRow(j,k) << ".000000000 ";
// 			else
// 				cout << thisCompressedSparseRow(j,k) << " ";
// 		}
// 		cout << "; ";
// 	}
// 	cout << "]";
// }


void CompressedSparseRow::show_vectors(){
	cout << "value:" << endl;
	for(int i =0; i<value_.size(); i++){
		cout << value_[i] << " " ;
	}
	cout << endl;
	cout << "row_ptr:" << endl;
		for(int i =0; i<row_ptr_.size(); i++){
		cout << row_ptr_[i] << " " ;
	}
	cout << endl;
	cout << "col_ind:" << endl;
		for(int i =0; i<col_ind_.size(); i++){
		cout << col_ind_[i] << " " ;
	}
	cout << endl;
}


void CompressedSparseRow::Show()
{
	//CompressedSparseRow& thisCSR = *this;
	if(value_.empty()){
		cout << "Matriz vacia" << endl;
	}else{
		for (int j = 0; j < rows_; j++) {
			if (row_ptr_[j] == -1){		//fila de ceros
				for(int i = 0; i < cols_; i++){
					cout << 0 << " ";
				}
			}else{
				int it_aux = row_ptr_[j]; //aca empieza la i-esima fila
				for(int i = 0; i < cols_; i++){
					if(i == col_ind_[it_aux]){
						cout << value_[it_aux] << " ";
						it_aux++;
					}else{
						cout << 0 << " ";
					}
				}
			}	
			cout << endl;
		}
	}
}



/*
void CompressedSparseRow::GetKeys(){
std::vector<char> v;
	for(pair<int,int>,double>::iterator it = data_.begin(); it != data_.end(); ++it) {
		v.push_back(it->first);
		cout << it->first << "\n";
	}
}
*/

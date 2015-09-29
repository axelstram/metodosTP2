#include "csr.h"


CompressedSparseRow::CompressedSparseRow(size_t n, size_t m) : rows_(n), cols_(m) {}


CompressedSparseRow::~CompressedSparseRow() {}


//usado para leer un elemento
double CompressedSparseRow::operator()(size_t i, size_t j)
{
	int k = row_ptr_[i];
	int condicion;
	if(i == rows_ -1){
		//soy la ultima fila, row_ptr_[i+1] se indefine 
		condicion = value_.size()-1;
	}else{
		condicion = row_ptr_[i+1]-1;
	}

	do{
		if(j == col_ind_[k]) return value_[k];
		k++;
	}while(k <= condicion);
	//no me encontre, soy cero
	return 0;
}


//para agregar un elemento
void CompressedSparseRow::operator()(size_t i, size_t j, double value)
{
	//matriz vacia
	if(value_.empty()){
		//si no es la primer fila
		if(i != 0){
			//filas vacias
			for(int k = 0; k < i; k++){
				value_.push_back(-1.0);
				col_ind_.push_back(-1);
				row_ptr_.push_back(k);
			}
			row_ptr_.push_back(i);
		}else{
			row_ptr_.push_back(0);
		}
		//termino de agregar el resto de la info
		value_.push_back(value);
		col_ind_.push_back(j);
	}else{

		if(i == row_ptr_.size()-1){
			//sigo agregando elementos a fila i
			col_ind_.push_back(j);
			value_.push_back(value);
		}else{
			if(i == row_ptr_.size()){
				//fila nueva
				col_ind_.push_back(j);
				value_.push_back(value);
				row_ptr_.push_back(value_.size()-1);
			}
			if(i > row_ptr_.size()){
				//agrego filas de ceros y luego la nueva fila
				for(int k = 0; k < i - row_ptr_.size(); k++){//agrego las filas ceros
					col_ind_.push_back(-1.0);
					value_.push_back(-1.0);
					row_ptr_.push_back(value_.size()-1);
				}
				//y agrego la nueva fila
				col_ind_.push_back(j);
				value_.push_back(value);
				row_ptr_.push_back(value_.size()-1);
			}
		}
	}
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
	CompressedSparseRow& thisCSR = *this;
	if(value_.empty()){
		cout << "Matriz vacia" << endl;
	}else{
		for (int i = 0; i < rows_; i++) {
			for(int j = 0; j < cols_; j++){
				cout << thisCSR(i, j) << " ";	
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

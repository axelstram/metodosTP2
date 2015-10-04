#include "csr.h"
#include <vector>
#include <cassert>

CompressedSparseRow::CompressedSparseRow(size_t n, size_t m) : rows_(n), cols_(m) {

	for(int i = 0; i<n; i++) row_ptr_.push_back(-1);
}

CompressedSparseRow::CompressedSparseRow(size_t n, size_t m, vector<double> value, vector<int> row_ptr, vector<int> col_ind)
{
	rows_ = n;
	cols_ = m;

	value_ = value; 
	row_ptr_ = row_ptr; 
	col_ind_ = col_ind; 
}

//usado para leer un elemento
double CompressedSparseRow::operator()(size_t i, size_t j) const
{
	int k = row_ptr_[i];
	int condicion;
	if(i < row_ptr_.size() -1){
		condicion = row_ptr_[i+1]-1;
	}else{
		if(i == row_ptr_.size()){
			//soy la ultima fila definida, row_ptr_[i+1] se indefine 
			condicion = value_.size()-1;
		}else{
			//caso en que tengo menos filas definidas que rows_(las filas posta)
			condicion = rows_;
		}
	}

	do{
		if(j == col_ind_[k]) return value_[k];
		k++;
	}while(k <= condicion);
	//no me encontre, soy cero
	return 0;
}


//para agregar/modificar un elemento, devuelve referencia
double& CompressedSparseRow::operator()(size_t i, size_t j)
{
	if(row_ptr_[i] == -1){
		//caso soy fila vacia, busco mi siguiente no vacio
		int next_row;
		for(next_row = i+1; next_row < rows_; next_row++){
			if(row_ptr_[next_row] != -1){
				break;
			}
		}
		if(next_row == rows_){
			//no habia ningun siguiente, soy el ultimo
			row_ptr_[i] = value_.size();
			col_ind_.push_back(j);
			value_.push_back(0);
			return value_.back();
		}else{
			//encontre la fila siguiente a mi
			int p = row_ptr_[next_row];
			value_.emplace(value_.begin() + p, 0);
			col_ind_.emplace(col_ind_.begin() + p, j);
			row_ptr_[i] = p;
			//corrijo el indice de mis siguientes validos
			for(int h = i+1; h < rows_; h++){
				if(row_ptr_[h] != -1) row_ptr_[h]++;
			}

			return value_.at(p);
		}
	}else{
		//caso mi fila no es vacia, me agrego en mi zona
		//me fijo si mi posicion es anterior al row_pointer
		if(j < col_ind_[row_ptr_[i]]){
			int p = row_ptr_[i];
			value_.emplace(value_.begin() + p, 0);
			col_ind_.emplace(col_ind_.begin() + p, j);

			//corrijo el indice de mis siguientes validos
			for(int h = i+1; h < rows_; h++){
				if(row_ptr_[h] != -1) row_ptr_[h]++;
			}

			return value_.at(p);
		}else{
			// busco de donde hasta donde son los values de mi row
			int from = row_ptr_[i];
			int next_row_ptr = value_.size(); 
			for (int h = i+1; h < row_ptr_.size(); ++h)
			{
				if(row_ptr_[h] != -1){
					next_row_ptr = row_ptr_[h];
					break;
				}
			}
			// from = indice del primer elem de la row en value
			// next_row_ptr = indice del primer elem de la sigueinte row no nula, a lo sumo es el final 
			for ( ; from <= next_row_ptr; ++from)
			{
				if(j == col_ind_[from]){
					return value_.at(from);
				}else{
					if(j < col_ind_[from] || from == next_row_ptr){
						// agrego antes de value_[from]
						value_.emplace(value_.begin() + from, 0);
						col_ind_.emplace(col_ind_.begin() + from , j);
						//corrijo el indice de mis siguientes validos
						for(int h = i+1; h < rows_; h++){
							if(row_ptr_[h] != -1) row_ptr_[h]++;
						}
						return value_.at(from);
					}
				}
					
			}
		}
	}
}


size_t CompressedSparseRow::rows() const
{
	return rows_;
}



size_t CompressedSparseRow::cols() const
{
	return cols_;
}

vector<double> CompressedSparseRow::values() const
{
	return value_;
}

vector<int> CompressedSparseRow::row_ptr() const
{
	return row_ptr_;
}

vector<int> CompressedSparseRow::col_ind() const
{
	return col_ind_;
}



CompressedSparseRow* CompressedSparseRow::operator*(double scalar)
{
	CompressedSparseRow& thisMat = *this;
	vector<double> newValues = thisMat.values();

	if (scalar != 0){
		for (int i = 0; i < newValues.size(); ++i){
			newValues[i] *= scalar;
		}
		CompressedSparseRow* res = new CompressedSparseRow(thisMat.cols(), thisMat.rows(), newValues, thisMat.row_ptr(), thisMat.col_ind());
		return res;

	}else{
		CompressedSparseRow* res = new CompressedSparseRow(thisMat.cols(), thisMat.rows());
		return res;
	}

}

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

vector<double> CompressedSparseRow::operator*(const vector<double>& x)
{
	vector<double> res(rows(), 0);

	for(int j = 0; j< rows(); j++){
		//si soy el ultimo la puta condicion se indefine
		int condicion;
		if(j == rows() -1){
			condicion = value_.size() -1;
		}else{
			condicion = row_ptr_[j+1]-1;
		}
	    for(int i = row_ptr_[j]; i <= condicion; i++){
	        res[j] += value_[i] * x[col_ind_[i]];
		}
	} 

	return res;  
}
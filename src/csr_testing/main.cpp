#include "csr.h"


int main(){

	CompressedSparseRow B = CompressedSparseRow(3, 3);
	B(0,0,1);
	B(1,1,1);
	B(2,2,1);
	B.Show();

	cout << endl; 

	CompressedSparseRow A = CompressedSparseRow(3, 3); 
	A(0,0,1);
	A(0,2,1);
	A(1,1,1);
	A(2,0,1);
	A(2,2,1);
	A.Show();

	return 0;
}
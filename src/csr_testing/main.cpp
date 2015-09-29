#include "csr.h"


int main(){

	CompressedSparseRow A = CompressedSparseRow(3, 3); 
	A(0,0) = 5.43;
	A(1,2) = 6;
	A(1,0) = 7;
	A(1,1) = 8;
	A(2,1) = 10;
	A(2,0) = 9;
	A(2,0) = 444;
	A(0,2) = 9.43;
	A(0,1) = 6.43;
	A.show_vectors();
	A.Show();
	//A.show_vectors();
	
	// cout << endl; 

	// CompressedSparseRow B = CompressedSparseRow(3, 3);
	// B(0,0) = 1;
	// B(1,1) = 1;
	// B(2,2) = 1;
	// B.Show();

	// cout << endl; 


	// CompressedSparseRow C = CompressedSparseRow(4, 4);
	// C(0,0) = 1;
	// C(1,1) = 1;
	// C(2,2) = 1;
	// C.Show();

	return 0;
}
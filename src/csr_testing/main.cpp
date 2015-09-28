#include "csr.h"


int main(){

	 CompressedSparseRow A = CompressedSparseRow(3, 3);
	 A.Show();
	 A(0,0,1);
	 A(1,1,1);
	 A(2,2,1);
	 A.Show();
	 A.show_vectors();

	return 0;
}
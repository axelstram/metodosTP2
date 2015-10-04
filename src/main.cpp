
#include "aux.h"

using namespace std;


///from test.in
int method;
double c; //prob. de teletransportacion
int instance_type;
string test_in_file;
string graph_file;
double tolerance;

///from test.txt
int nodes; /* 	 pages / teams	 	*/
int edges; /*	 links / marches 	*/
int matrix_type = DOK_MATRIX;


int main(int argc, char* argv[])
{	
 	Matrix& A = load_test_in(argv[1]);

//-----------test-------------
/*

	///generar c aleatorios
	UniformDist udist(0,1);

	vector<double> cs;
	for(int i=0;i<100;i++){
		double c =udist.Generate();

		//truncar a 2 decimales
		c *= 100.;
		int ci = c;
		c = (double)ci/100.;

		cs.push_back(c);
		cout<<c<<endl;
	}
*/
//-----------------------

	vector<double> x(A.cols());
	for (int i = 0; i < x.size(); i++)
		x[i] = 1;

	int maxIter = 200000;
	pair<double, vector<double>> res;
	pair<double, vector<double>> res2;

	bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

	if (encontroResultado) {
		show_vector(res.second);
		escribir_resultado(res.second, argv[2]);
	} else {
 		cout << "no encontro resultado" << endl;
    }

    delete &A;
}




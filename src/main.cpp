
#include "aux.h"
#include <chrono>


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

bool medir_tiempo = false;

int main(int argc, char* argv[])
{	
	if(argc > 3)
		medir_tiempo = true;


	vector<long int> times;
	for (int t = 0; t < 10; ++t)
	{
	 	Matrix& A = load_test_in(argv[1]);

		vector<double> x(A.cols());
		for (int i = 0; i < x.size(); i++)
			x[i] = 1;

		int maxIter = 200000;
		pair<double, vector<double>> res;
		pair<double, vector<double>> res2;

		long int total_time = 0;
		auto begin = std::chrono::high_resolution_clock::now();

		bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

		auto end = std::chrono::high_resolution_clock::now();

		auto power_method_time = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
		total_time += power_method_time ; 
		times.push_back(total_time);
	
		if (encontroResultado) {
			//show_vector(res.second);
			escribir_resultado(res.second, argv[2]);
		} else {
	 		cout << "no encontro resultado" << endl;
	    }
	
  		
  		delete &A;

  		if (!medir_tiempo)
  			return 0;
	}

	long int promedio = 0; 
	for (int i = 0; i < times.size(); ++i)
	{	
		promedio += times[i];
	}

	cout << promedio/times.size() << endl; // en ms

	return 0;
}




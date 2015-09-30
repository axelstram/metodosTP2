#include <string>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Mat.h"

#define PAGE_RANK_METHOD 0
#define ALT_METHOD 1

#define WEB_RANK 0
#define SPORT_RANK 1

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

void show_vector(vector<double> v){
	for(int i = 0; i < v.size(); i++){
		cout << v[i] << ' ';
	}
	cout << endl;
} 



Mat LinkMatrixModification(Mat A, double m)
{
	Mat S(A.cols(), A.rows());

	for (int i = 0; i < S.rows(); i++) {
		for (int j = 0; j < S.cols(); j++) {
			S(i, j) = 1.0/(double)S.rows();
		}
	}

	Mat M(A.cols(), A.rows());

	M = A*(1.0-m) + S*m;

	return M;
}




Mat loadWebGraph(string graph_file)
{
	ifstream f(graph_file);
	string s;
	//int nodes;
	//int edges;

	//ignoramos todo hasta llegar a nodes
	while (s != "Nodes:")
		f >> s;

	f >> s;
	nodes = stoi(s);

	f >> s; //Edges:
	f >> s;
	edges = stoi(s);

	//ignoro todo hasta ToNodeId
	while (s != "ToNodeId")
		f >> s;

	vector<int> linksEntrantes(nodes);
	Mat A(nodes, nodes);

	for (int i = 0; i < edges; i++) {
		int from;
		int to;

		f >> s;
		from = stoi(s);
		f >> s;
		to = stoi(s);
		from--;
		to--;

		linksEntrantes[from]++;
		A(to, from) = 1;
	}

	for (int from = 0; from < A.rows(); from++) {
		for (int to = 0; to < A.cols(); to++) {
			A(to, from) = A(to, from)/(double)linksEntrantes[from];
		}
	}


	return A;
}


void normalizarMatrizEquipos( Mat& A){

	//para cada elemento i j....
	for (int i = 0; i < A.rows(); ++i){
		//sumo la columna i
		double acum = 0;
		for (int k = 0; k < A.cols(); ++k){
			acum += A(k,i);
		}
		//divido cada elemento de la fila i por acum
		for (int j = 0; j < A.cols(); ++j){
			A(i,j) /= acum;
		}
	}
}


Mat loadSportGraph(string graph_file){
	Mat A(nodes,nodes);

	ifstream f(graph_file);
	string s;

	f >> s;
	nodes = stoi(s);

	f >> s;
	edges = stoi(s);

	for (int i = 0; i < edges; i++) {
		int fecha;
		int equipo1,equipo2,goles_equipo1,goles_equipo2;

		f >> s;
		fecha = stoi(s);
		f >> s;
		equipo1 = stoi(s);
		f >> s;
		goles_equipo1 = stoi(s);
		f >> s;
		equipo2 = stoi(s);
		f >> s;
		goles_equipo2 = stoi(s);

		equipo1--;
		equipo2--;


		if(goles_equipo2>goles_equipo1)
			A(equipo1, equipo2) += goles_equipo2-goles_equipo1;

		if(goles_equipo1>goles_equipo2)
			A(equipo2, equipo1) += goles_equipo1-goles_equipo2;

	}

	///normalizar matriz..
	normalizarMatrizEquipos(A);



	return A;
}


Mat load_test_in(string test_in_file){

	ifstream f(test_in_file);
	string s;
	

	f >> s;
	method = stoi(s);
	f >> s;
	c = stod(s);
	f >> s;
	instance_type = stoi(s);
	f >> s;
	graph_file = s;
	f >> s;
	tolerance = stod(s);

	///load matrix;
	if(instance_type == WEB_RANK){
		Mat A = loadWebGraph(graph_file);
		return A;
	}else{
		if(instance_type == SPORT_RANK){
			cout << "aca" << endl;
			Mat A = loadSportGraph(graph_file);
		}
	}
}


// Nuestro load
Mat LoadMatrixFromFile(string file_path)
{
	ifstream f(file_path);
	string s;
	
	f >> s;
	int dim = stoi(s);
	
	Mat A(dim, dim);

	if (dim != A.rows()) {
		cout << "Dimensiones incorrectas al cargar " << endl << "dim " << dim << endl << "rows(A) " << A.rows() << endl;
		exit(1);
	}

	//Cargo A
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			f >> s;
			A(i,j) = stod(s);
		}
	}

	return A;
}



double normaInfVec(vector<double> v){
	double res = 0;
	int indice = 0;
	
	for (int i = 0; i < v.size(); i++) {
		if (fabs(v[i]) > fabs(res)) 
			res = fabs(v[i]);
	}

	return res;
}



vector<double> vec_sub(vector<double>& x, vector<double>& y)
{
	vector<double> resta(x.size());

	assert(x.size() == y.size());

	for (int i = 0; i < x.size(); i++)
		resta[i] = x[i] - y[i];

	return resta;
}



bool MetodoPotencia(Mat& A, vector<double> x,double c, float tolerance, int maxIter, pair<double, vector<double>>& res)
{
	int k = 1;
	double infNormX = normaInfVec(x);
	//double lambda = 0;
	double anterior;

		double z = 0;
		for(int i = 0;i<nodes;i++) z += c*x[i]*1./(double)nodes;

		for(int i = 0;i<nodes;i++) x[i] = x[i]*(1.-c)+z;

	
	for (int i = 0; i < x.size(); i++) {
		x[i] /= infNormX;
	}
	

	while (k <= maxIter) {
		anterior = infNormX;
		
		vector<double> y = A*x;
		//vec_sum

		z=0;
		for(int i = 0;i<nodes;i++) z += c*y[i]*1./(double)nodes;

		for(int i = 0;i<nodes;i++) y[i] = y[i]*(1.-c) + z;

		double infNormY = normaInfVec(y);		
	
		
		if (infNormY == 0) {
			cout << "Vector inicial incorrecto" << endl;
			exit(1);
		}
		
		for (int i = 0; i < y.size(); i++) {
			y[i] /= infNormY;
		}
		
		//show_vector(y);


		
		double error = normaInfVec(vec_sub(x, y));
		
		//double error = fabs(lambda - para);
		
		if (error < tolerance) {
			res = make_pair(infNormY, y);
			return true;
		}

		x = y;

		k++;
	}

	return false;
}

bool MetodoPotencia2(Mat& A, vector<double> x, float tolerance, int maxIter, pair<double, vector<double>>& res);

int main(int argc, char* argv[])
{


	
	//Mat A = LoadMatrixFromFile("/home/sebs/Desktop/metodosTP2/src/multPrueba.txt");	

	Mat A = load_test_in(argv[1]);
	A.Show();
	/*Mat M = LinkMatrixModification(A, c);

	vector<double> x = {1, 1, 1, 1};
	int maxIter = 10000000;
	pair<double, vector<double>> res;
	pair<double, vector<double>> res2;

	bool encontroResultado2 = MetodoPotencia2(M, x, tolerance, maxIter, res2);

	bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);



	if (encontroResultado) {
		cout << "autovalor " << res.first << endl; 
		cout << "autovector " << endl;
		show_vector(res.second);
		cout << endl;

		vector<double> Mx = M*res.second;
		vector<double> lambda_x(Mx.size());

		for (int i = 0; i < lambda_x.size(); i++)
			lambda_x[i] = res.first * res.second[i];

		vector<double> Mx_menos_lambda_x = vec_sub(Mx, lambda_x);

		for (int i = 0; i < Mx_menos_lambda_x.size(); i++) {
			if (Mx_menos_lambda_x[i] > 0.001) {
				cout << "Mx y lambda*x distintos" << endl;
				cout << "Mx: ";
				show_vector(Mx);
				cout << "lambda*x: ";
				show_vector(lambda_x);
			
				return 0;
			}
		}

		cout << "Mx y lambda*x iguales" << endl;		
	} else {
		cout << "no encontro resultado" << endl;
	}
	*/
}







bool MetodoPotencia2(Mat& A, vector<double> x, float tolerance, int maxIter, pair<double, vector<double>>& res)
{
	int k = 1;
	double infNormX = normaInfVec(x);
	//double lambda = 0;
	double anterior;

	
	for (int i = 0; i < x.size(); i++) {
		x[i] /= infNormX;
	}
	

	while (k <= maxIter) {
		anterior = infNormX;
		vector<double> y = A*x;

		double infNormY = normaInfVec(y);		
	
		
		if (infNormY == 0) {
			cout << "Vector inicial incorrecto" << endl;
			exit(1);
		}
		
		for (int i = 0; i < y.size(); i++) {
			y[i] /= infNormY;
		}
		
		//show_vector(y);


		
		double error = normaInfVec(vec_sub(x, y));
		
		//double error = fabs(lambda - para);
		
		if (error < tolerance) {
			res = make_pair(infNormY, y);
			return true;
		}

		x = y;

		k++;
	}

	return false;
}
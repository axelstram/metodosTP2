#include <string>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "Mat.h"
#include "dok.h"
#include "csr.h"


#define PAGE_RANK_METHOD 0
#define ALT_METHOD 1

#define WEB_RANK 0
#define SPORT_RANK 1


#define VECTOR_MATRIX 0
#define DOK_MATRIX 1
#define CSR_MATRIX 2

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
int matrix_type = CSR_MATRIX;

void show_vector(vector<double> v){
	for(int i = 0; i < v.size(); i++){
		cout << v[i] << ' ';
	}
	cout << endl;
} 




vector<double> vec_sum(const vector<double>& x, const vector<double> y){
	vector<double> res;
	for(unsigned int i = 0; i < x.size(); i++){
		res.push_back(x[i] + y[i]);
	}
	return res;
}



vector<double> vec_mult_scalar(const vector<double>& x, const double scalar)
{
	vector<double> res;
	for(unsigned int i = 0; i < x.size(); i++){
		res.push_back(x[i] * scalar);
	}
	return res;	
}

/*
Matrix LinkMatrixModification(Matrix A, double m)
{
	Matrix S(A.cols(), A.rows());

	for (int i = 0; i < S.rows(); i++) {
		for (int j = 0; j < S.cols(); j++) {
			S(i, j) = 1.0/(double)S.rows();
		}
	}

	Matrix M(A.cols(), A.rows());

	M = A*(1.0-m) + S*m;

	return M;
}
*/



Matrix& loadWebGraph(string graph_file)
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

	Matrix* A_tmp;
	switch ( matrix_type )
    {
        case VECTOR_MATRIX :
            A_tmp = new Mat(nodes, nodes);
            break;
        case DOK_MATRIX :
            A_tmp = new DictionaryOfKeys(nodes, nodes);
            break;
        case CSR_MATRIX :
         	A_tmp = new CompressedSparseRow(nodes,nodes);
         	break;
    }
    Matrix& A = *A_tmp;



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


void normalizarMatrizEquipos( Matrix& A){

	//para cada elemento i j....
	for (int i = 0; i < A.cols(); ++i){
		//sumo la columna i
		double acum = 0;
		for (int k = 0; k < A.rows(); ++k){
			acum += A(k,i);
		}

		//divido cada elemento de la fila i por acum
		for (int j = 0; j < A.rows(); ++j){
			if (acum != 0)
				A(j,i) /= acum;
			else
				A(j,i) = 1./A.rows();
		}
	}
}


Matrix& loadSportGraph(string graph_file){


	ifstream f(graph_file);
	string s;

	f >> s;
	nodes = stoi(s);

	f >> s;
	edges = stoi(s);


//	Matrix A(nodes,nodes);
	Matrix* A_tmp;

	switch ( matrix_type )
    {
        case VECTOR_MATRIX :
            A_tmp = new Mat(nodes, nodes);
            break;
        case DOK_MATRIX :
            A_tmp = new DictionaryOfKeys(nodes, nodes);
            break;
        case CSR_MATRIX :
         	A_tmp = new CompressedSparseRow(nodes,nodes);
         	break;
    }
    Matrix& A = *A_tmp;

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
			A(equipo2, equipo1) += (double)(goles_equipo2-goles_equipo1);

		if(goles_equipo1>goles_equipo2)
			A(equipo1, equipo2) += (double)(goles_equipo1-goles_equipo2);

	}
	A.Show();
	cout<<endl<<endl;
	///normalizar matriz..
	normalizarMatrizEquipos(A);

	return A;
}


Matrix& load_test_in(string test_in_file){

	ifstream f(test_in_file);
	string s;
	

	f >> s;
	method = stoi(s);
	f >> s;
	c = stod(s);
	c = 1-c; //lo doy vuelta para que sea CONSISTENTE CON EL PAPER.
	f >> s;
	instance_type = stoi(s);
	f >> s;
	graph_file = s;
	f >> s;
	tolerance = stod(s);

	///load matrix;
	if(instance_type == WEB_RANK){
		Matrix& A = loadWebGraph(graph_file);
		return A;

	}else if(instance_type == SPORT_RANK) {
		Matrix& A = loadSportGraph(graph_file);
		return A;
	}

}


// Nuestro load
Matrix& LoadMatrixFromFile(string file_path)
{
	ifstream f(file_path);
	string s;
	
	f >> s;
	int dim = stoi(s);
	
	//Matrix A(dim, dim);
	Matrix* A_tmp;

	switch ( matrix_type )
    {
        case VECTOR_MATRIX :
            A_tmp = new Mat(nodes, nodes);
            break;
        case DOK_MATRIX :
            A_tmp = new DictionaryOfKeys(nodes, nodes);
            break;
        case CSR_MATRIX :
         	A_tmp = new CompressedSparseRow(nodes,nodes);
         	break;
    }
    Matrix& A = *A_tmp;

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


double normaUnoVec(vector<double> v){
	double res = 0;
	
	for (int i = 0; i < v.size(); i++) {
			res += fabs(v[i]);
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



bool MetodoPotencia(Matrix& A, vector<double> x,double c, float tolerance, int maxIter, pair<double, vector<double>>& res)
{
	int k = 1;
	double NormX = normaUnoVec(x);

	for (int i = 0; i < x.size(); i++) {
		x[i] /= NormX;
	}

	//double ms = c/(double)nodes;
	vector<double> ms(A.rows());
	for (int i = 0; i < ms.size(); i++)
		ms[i] = c/(double)nodes;

	while (k <= maxIter) {		
		//M = A*(1.0-m) + m*S; 
		/*
		Matrix A2 = A*(1.0-c);
		//sumo a todas las posiciones de A2 el escalar ms.
		//es mas eficiente en cuanto a espacio que crear la matriz m*S explicitamente y sumarsela a A2.
		Matrix A2_mas_ms = A2 + ms;

		vector<double> y = A2_mas_ms*x;		
		*/
		Matrix* A_prima = A*(1.0-c);

		vector<double> y = vec_sum((*A_prima)*x, ms);
		delete(A_prima);

		double normY = normaUnoVec(y);			
		
		if (normY == 0) {
			cout << "Vector inicial incorrecto" << endl;
			exit(1);
		}

		for (int i = 0; i < y.size(); i++) {
			y[i] /= normY;
		}
		
		double error = normaUnoVec(vec_sub(x, y));

		if (error < tolerance) {
			res = make_pair(normY, y);
			return true;
		}

		x = y;

		k++;
	}

	return false;
}


/*

void ChequearAutovector(Mat& A, double autoval, vector<double>& x)
{
	cout << "autovalor " << autoval << endl; 
	cout << "autovector " << endl;
	show_vector(x);
	cout << endl;

	vector<double> Ax = A*x;
	vector<double> lambda_x(Ax.size());

	for (int i = 0; i < lambda_x.size(); i++)
		lambda_x[i] = autoval * x[i];

	vector<double> Ax_menos_lambda_x = vec_sub(Ax, lambda_x);

	for (int i = 0; i < Ax_menos_lambda_x.size(); i++) {
		if (Ax_menos_lambda_x[i] > 0.001) {
			cout << "Ax y lambda*x distintos" << endl;
			cout << "Ax: ";
			show_vector(Ax);
			cout << "lambda*x: ";
			show_vector(lambda_x);
		
			return;
		}
	}

	cout << "Ax y lambda*x iguales" << endl;
}

*/

bool comparePair(pair<int,int> p1,pair<int,int> p2)
{
	return p1.second > p2.second;
}



vector<pair<int,int> > IN_DEG(Matrix& A)
{

	vector<pair<int,int> > rank;

	for(int i = 0;i<A.rows();i++){
		int links_entrantes = 0;
		for(int j = 0;j<A.cols();j++){
			if(A(i,j)!=0)links_entrantes++;
		}
		rank.push_back(make_pair(i,links_entrantes));
	}

	sort(rank.begin(), rank.end(), comparePair);

	for(int i=0;i<rank.size();i++)
		cout<< "("<<rank[i].second<<","<<rank[i].first<<") ";
	cout<<endl;

	return rank;

}



void escribir_resultado(vector<double>& x, string output_path)
{
	ofstream outputFile(output_path);

	for (int i = 0; i < x.size(); i++)
		outputFile << x[i] << endl; 
}


//bool MetodoPotencia2(Mat& A, vector<double> x, float tolerance, int maxIter, pair<double, vector<double>>& res);
//vector<double> power_method_d(Matrix A,vector<double> v);
//void power_method_short(Matrix A,vector<double> x);

int main(int argc, char* argv[])
{	
	//Matrix A = LoadMatrixFromFile("/home/vivi/metodosTP2/src/matrizpiolaM.txt");	
 	Matrix& A = load_test_in(argv[1]);
 	//Matrix M = LinkMatrixModification(A, c);

	//IN_DEG(A);

	vector<double> x(A.cols());
	for (int i = 0; i < x.size(); i++)
		x[i] = 1;

	int maxIter = 200000;
	pair<double, vector<double>> res;
	pair<double, vector<double>> res2;

	//bool encontroResultado = MetodoPotencia2(M, x, tolerance, maxIter, res);
 
	bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

	if (encontroResultado) {
		show_vector(res.second);
		escribir_resultado(res.second, argv[2]);
	} else {
 		cout << "no encontro resultado" << endl;
    }

    delete &A;
}




/*


bool MetodoPotencia2(Mat& A, vector<double> x, float tolerance, int maxIter, pair<double, vector<double>>& res)
{
	int k = 1;
	double NormX = normaUnoVec(x);
	
	for (int i = 0; i < x.size(); i++) {
		x[i] /= NormX;
	}
	

	while (k <= maxIter) {
		vector<double> y = A*x;

		double NormY = normaUnoVec(y);		
	
		if (NormY == 0) {
			cout << "Vector inicial incorrecto" << endl;
			exit(1);
		}
		
		for (int i = 0; i < y.size(); i++) {
			y[i] /= NormY;
		}	
		
		double error = normaUnoVec(vec_sub(x, y));
			
		if (error < tolerance) {
			res = make_pair(NormY, y);
			return true;
		}

		x = y;

		k++;
	}

	return false;
}



*/
#include <string>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <algorithm>
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


Mat loadSportGraph(string graph_file){


	ifstream f(graph_file);
	string s;

	f >> s;
	nodes = stoi(s);

	f >> s;
	edges = stoi(s);


	Mat A(nodes,nodes);


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

	}else if(instance_type == SPORT_RANK) {
		Mat A = loadSportGraph(graph_file);
		return A;
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
	//double lambda = 0;
	double anterior;


	
	double infNormX = normaInfVec(x);

	for (int i = 0; i < x.size(); i++) {
		x[i] /= infNormX;
	}
	
	double z = 0;
	for(int i = 0;i<nodes;i++) z += c*x[i]*1./(double)nodes;

	for(int i = 0;i<nodes;i++) x[i] = x[i]*(1.-c);




	while (k <= maxIter) {
		anterior = infNormX;
		
		vector<double> y = A*x;
		for (int i = 0; i < y.size(); i++)
			y[i] += z;


		double infNormY = normaInfVec(y);			
		
		if (infNormY == 0) {
			cout << "Vector inicial incorrecto" << endl;
			exit(1);
		}
		
		for (int i = 0; i < y.size(); i++) {
			y[i] /= infNormY;
		}
		


		z = 0;
		for(int i = 0;i<nodes;i++) z += c*y[i]*1./(double)nodes;

		for(int i = 0;i<nodes;i++) y[i] = y[i]*(1.-c);

		//show_vector(y);


		
		double error = normaInfVec(vec_sub(x, y));
				
		if (error < tolerance) {
			res = make_pair(infNormY, y);
			return true;
		}

		x = y;

		k++;
	}

	return false;
}




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

bool comparePair(pair<int,int> p1,pair<int,int> p2){

return p1.second > p2.second;
}

vector<pair<int,int> > IN_DEG(Mat A)
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

for(int i=0;i<rank.size();i++)cout<< "("<<rank[i].second<<","<<rank[i].first<<") ";
	cout<<endl;

	return rank;

}





bool MetodoPotencia2(Mat& A, vector<double> x, float tolerance, int maxIter, pair<double, vector<double>>& res);
vector<double> power_method_damian(Mat A,vector<double> v);
void power_method_internet_v1(Mat A,vector<double> x);

int main(int argc, char* argv[])
{	



	//Mat A = LoadMatrixFromFile("/home/vivi/metodosTP2/src/matrizpiolaM.txt");	
 	Mat A = load_test_in(argv[1]);
 	Mat M = LinkMatrixModification(A, c);

cout<<"------------------"<<endl;


IN_DEG(A);


cout<<"------------------"<<endl;
 	A.Show();
 	cout<<endl;
 	M.Show();



 	vector<double> x = {1, 1, 1, 1};
power_method_internet_v1(A,x);

/*

	//A.Show();
	
	Mat M = LinkMatrixModification(A, c);
	M.Show();
	vector<double> x = {1, 1, 1, 1};
	int maxIter = 10000000;
	pair<double, vector<double>> res;
	pair<double, vector<double>> res2;

	bool encontroResultado = MetodoPotencia2(M, x, tolerance, maxIter, res2);

	//bool encontroResultado2 = MetodoPotencia(A, x, c , tolerance, maxIter, res2);

//res2.second = power_method_damian(A,x);

		if (encontroResultado) {

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
		//x[i] /= infNormX;
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








///////////////////////////////////////////

///////////////////////////////////////////

vector<double> restaVectores(vector<double>& y,vector<double>& x){
	vector<double> res;
	for(unsigned int i = 0; i < y.size(); i++){
		res.push_back(y[i] - x[i]);
	}
	return res;
}


vector<double> vectorXescalar(vector<double>& v , double w){
	vector<double> res;
	for(unsigned int i = 0; i < v.size(); i++){
		res.push_back(v[i] * w);
	}
	return res;	
}


vector<double> sumaVectores(vector<double>& x,vector<double> y){
	vector<double> res;
	for(unsigned int i = 0; i < y.size(); i++){
		res.push_back(x[i] + y[i]);
	}
	return res;
}


double norm_uno(vector<double> y){
	double res = 0;
	for(unsigned int i = 0; i < y.size(); i++){
		res = res + fabs(y[i]);
	}
	return res;
}

vector<double> power_method_damian(Mat A,vector<double> v){
	double delta;
	vector<double> x = v;
	vector<double> y;
	double w;

	int i = 0;
	do{
		y = A*x;
		w = norm_uno(x) - norm_uno(y);//norma inf
		y = sumaVectores(y,vectorXescalar(v,w)); 
		
		delta = norm_uno(restaVectores(y,x));
		x = y;
		i++;
	}while(!(delta < tolerance));
	//~ cout << i << endl;
	return x;
}













void power_method_internet_v1(Mat A,vector<double> x)
{
double temp;
int n = nodes;
int d=0;
vector<double> cv;
do
    {
        cv = A*x;
        for(int i=0;i<n;i++)
            x[i]=cv[i];
            
        temp=d;
        d=0;
        
        for(int i=0;i<n;i++)
        {
            if(fabs(x[i])>fabs(d))
                d=x[i];
        }

        for(int i=0;i<n;i++)
            x[i]/=d;
            
    }while(fabs(d-temp)>0.00001);


cout << "autovalor "<<d<<endl;

cout << "autovector "<<endl;
  for(int i=0;i<n;i++)
        cout<<" "<<x[i];

    cout<<endl;

    ChequearAutovector(A, d, x);
}
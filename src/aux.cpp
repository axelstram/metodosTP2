#include "aux.h"

extern int method;
extern double c; //prob. de teletransportacion
extern int instance_type;
extern string test_in_file;
extern string graph_file;
extern double tolerance;

///from test.txt
extern int nodes; /* 	 pages / teams	 	*/
extern int edges; /*	 links / marches 	*/
extern int matrix_type;


// funciones de vectores

// muestra un vector por pantalla
void show_vector(vector<double> v){
	for(int i = 0; i < v.size(); i++){
		cout << v[i] << ' ';
	}
	cout << endl;
} 

// suma de vectores
vector<double> vec_sum(const vector<double>& x, const vector<double> y){
	vector<double> res;
	for(unsigned int i = 0; i < x.size(); i++){
		res.push_back(x[i] + y[i]);
	}
	return res;
}

// norma uno de vector
double normaUnoVec(vector<double> v){
	double res = 0;
	
	for (int i = 0; i < v.size(); i++) {
			res += fabs(v[i]);
	}

	return res;
}

// resta de vectores
vector<double> vec_sub(vector<double>& x, vector<double>& y)
{
	vector<double> resta(x.size());

	assert(x.size() == y.size());

	for (int i = 0; i < x.size(); i++)
		resta[i] = x[i] - y[i];

	return resta;
}

// Funcion que normaliza lso elemntos de una matriz
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


// Funciones de cargado de matriz:

// Carga la matriz para paginas web
Matrix& loadWebGraph(string graph_file)
{
	ifstream f(graph_file);
	string s;

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

    cout << edges << endl;
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
			if(linksEntrantes[from] != 0)
				A(to, from) = A(to, from)/(double)linksEntrantes[from];
		}
		cout << from << endl;
	}
	return A;
}



// Carga la matriz para equipos
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
	//A.Show();

	cout<<endl<<endl;
	///normalizar matriz..
	normalizarMatrizEquipos(A);

	return A;
}


// Carga la matriz llamando al cargador de paginas web o al de deportes
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
		cout << "termino de cargar" << endl;
		return A;

	}else if(instance_type == SPORT_RANK) {
		Matrix& A = loadSportGraph(graph_file);
		return A;
	}

}

/*
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
*/






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

		// lo hacemos por paso porque devuelve un puntero a Matrix
		// Matrix es clase abstracta y no deja devolver por copia
		// entonces hacemos new y tenemos que hacer delete
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

		//double errorent=error*10000.;
		//cout<< (int) errorent<<endl;


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

// Comparador de Pairs para el sort del IN-DEG
bool comparePair(pair<int,int> p1,pair<int,int> p2)
{
	return p1.second > p2.second;
}

// Algoritmo IN-DEG
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

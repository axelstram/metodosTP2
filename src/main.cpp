
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



void ProcesarNormalmente(string input_file, string output_file)
{
 	Matrix& A = load_test_in(input_file);

	vector<double> x(A.cols(), 1);

	int maxIter = 200000;
	pair<double, vector<double>> res;

	bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

	if (encontroResultado) {
		escribir_resultado(res.second, output_file);
	} else {
 		cout << "no encontro resultado" << endl;
    }
		
	delete &A;
}



void MedirTiempos(string input_file, string output_file)
{
	vector<long int> times;
	Matrix& A = load_test_in(input_file);

	for (int t = 0; t < 10; ++t)
	{
		vector<double> x(A.cols(), 1);

		int maxIter = 200000;
		pair<double, vector<double>> res;

		long int total_time = 0;
		auto begin = std::chrono::high_resolution_clock::now();

		bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

		auto end = std::chrono::high_resolution_clock::now();

		auto power_method_time = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
		total_time += power_method_time ; 
		times.push_back(total_time);
	
		if (encontroResultado) {
			escribir_resultado(res.second, output_file);
		} else {
	 		cout << "no encontro resultado" << endl;
	    }
	}

	delete &A;

	long int promedio = 0; 
	sort(times.begin(), times.end());

	cout << times[times.size()/2] << endl;
}



//En input_file tiene que figurar el path a la carpeta que contiene todas las instancias que hay que correr.
void ProcesarBatchYMedirTiempos(string input_file, string output_file)
{
	ifstream inputFile(input_file);
	string s;

	inputFile >> s;
	method = stoi(s);
	inputFile >> s;
	c = stod(s);
	c = 1-c; //lo doy vuelta para que sea CONSISTENTE CON EL PAPER.
	inputFile >> s;
	instance_type = stoi(s);
	//Tomo el path
	inputFile >> s;
	string input_dir = s;
	inputFile >> s;
	tolerance = stod(s);

	inputFile.close();

	//Asumo que estoy parado en src
	string command = "ls -v " + input_dir + " > batch_file_names.txt";
	system(command.c_str());

	std::ifstream infile("batch_file_names.txt");
    string path;
	vector<string> batch;

	//guardo los paths de cada archivo individual.
    while(infile >> path){
        batch.push_back(input_dir + "/" + path);
    }

    for (int j = 0; j < batch.size(); j++) {
		vector<long int> times;

		Matrix& A = load_test_in_batch(batch[j]);

		for (int t = 0; t < 10; t++) {
			vector<double> x(A.cols(), 1);

			int maxIter = 200000;
			pair<double, vector<double>> res;

			long int total_time = 0;
			auto begin = std::chrono::high_resolution_clock::now();

			bool encontroResultado = MetodoPotencia(A, x, c , tolerance, maxIter, res);

			auto end = std::chrono::high_resolution_clock::now();

			auto power_method_time = chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
			total_time += power_method_time ; 
			times.push_back(total_time);

			if (encontroResultado) {
				escribir_resultado(res.second, output_file);
			} else {
		 		cout << "no encontro resultado" << endl;
		    }
		}

		delete &A;

		long int promedio = 0; 
		sort(times.begin(), times.end());

		cout << j+10 << " " << times[0] << endl;
	}

}



//1er parametro: archivo de entrada
//2do parametro: archivo de salida
//3er parametro: medir o no el tiempo
//4to parametro: procesar un batch
int main(int argc, char* argv[])
{	
	if (argc == 3)
		ProcesarNormalmente(argv[1], argv[2]);
	
	if (argc == 4)
		MedirTiempos(argv[1], argv[2]);

	if (argc >= 5)
		ProcesarBatchYMedirTiempos(argv[1], argv[2]);
	/*
	if (argc > 6)
		cout << "Cantidad de parametros incorrecta" << endl; exit(1);
*/

	return 0;
}




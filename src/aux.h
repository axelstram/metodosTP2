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
#include "rand.h"

#define PAGE_RANK_METHOD 0
#define ALT_METHOD 1

#define WEB_RANK 0
#define SPORT_RANK 1

// que tipo de matriz es la subyacente
#define VECTOR_MATRIX 0
#define DOK_MATRIX 1
#define CSR_MATRIX 2

using namespace std;


void show_vector(vector<double> v);
vector<double> vec_sum(const vector<double>& x, const vector<double> y);
double normaUnoVec(vector<double> v);
vector<double> vec_sub(vector<double>& x, vector<double>& y);
Matrix& loadWebGraph(string graph_file);
Matrix& loadSportGraph(string graph_file);
Matrix& load_test_in(string test_in_file);
Matrix& load_test_in_batch(string batch_instance_file);
void normalizarMatrizEquipos( Matrix& A);
bool MetodoPotencia(Matrix& A, vector<double> x,double c, float tolerance, int maxIter, pair<double, vector<double>>& res);
bool comparePair(pair<int,int> p1,pair<int,int> p2);
vector<pair<int,int> > IN_DEG(Matrix& A);
void escribir_resultado(vector<double>& x, string output_path);
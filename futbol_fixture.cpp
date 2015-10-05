#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

bool comparePair(pair<int,int> p1,pair<int,int> p2){
	return p1.second > p2.second;
}

void fixtureFutbol(string graph_file){

	ifstream f(graph_file);
	string s;
	string fecha;
	int eq1;
	int eq2;
	int goles_eq1;
	int goles_eq2;
	f >> s;
	int num_equipos = stoi(s);
	f >> s;
	int num_partidos = stoi(s);
	//cout << num_equipos << " " << num_partidos << endl;

	vector<pair<int,int> > equipo_puntaje (num_equipos+1);
	for (int i = 1; i < num_equipos+1; ++i){
		equipo_puntaje[i] = make_pair(i,0); //equipo i, 0 puntos inicialmente
	}

	for (int i = 0; i < num_partidos; ++i)
	{
			f >> fecha;
			f >> s;
			eq1 = stoi(s);
			
			f >> s;
			goles_eq1 = stoi(s);
			
			f >> s;
			eq2 = stoi(s);

			f >> s;
			goles_eq2 = stoi(s);
			
//			cout << fecha << " " << eq1 << " " << goles_eq1 << " " << eq2 << " " << goles_eq2 << endl;
			if(goles_eq1 < goles_eq2)
				equipo_puntaje[eq2].second += 3;
			
			if (goles_eq1 > goles_eq2)
				equipo_puntaje[eq1].second += 3;
			
			if (goles_eq1 == goles_eq2){
				equipo_puntaje[eq1].second += 1;
				equipo_puntaje[eq2].second += 1;
			}

	}

	sort(equipo_puntaje.begin()+1, equipo_puntaje.end(), comparePair);

	for (int i = 1; i < equipo_puntaje.size(); ++i){
		cout << equipo_puntaje[i].first << " " << equipo_puntaje[i].second << endl;
	}
}


int main(int argc, char* argv[]){
	fixtureFutbol(argv[1]);
}
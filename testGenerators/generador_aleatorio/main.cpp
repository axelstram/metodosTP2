#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "rand.h"

using namespace std;


bool isPair(pair<int,int> edge, vector<pair<int,int> > edges){
	for(int i=0;i<edges.size();i++){
		if(edge.first==edges[i].first &&edge.second==edges[i].second)
			return true;
	}
	return false;
}



void show_pair(pair<int,int> pair){
	cout<<pair.first<<" "<<pair.second;
}

void show_pair_vector(vector<pair<int,int> > v){
	for(int i=0;i<v.size();i++){
		show_pair(v[i]);
		cout<<endl;
	}
}



void getNodesAndEdges(string input_file, int& nodes, int& edges)
{
	ifstream inputStream(input_file);	
	string s;

	while (inputStream >> s) {
		nodes++;

		do {
			inputStream >> s;
			if (s == "-1")
				continue;
			edges++;
		} while (s != "-1");
	}

	inputStream.close();
}


void modifyInputGraph(string input_file, string output_file){
		int nodes = 0;
		int edges = 0;

		getNodesAndEdges(input_file, nodes, edges);

		ifstream inputStream(input_file);
		ofstream outputStream(output_file);
		string s;

		outputStream << nodes << " " << edges << endl;

		while (inputStream >> s) {
			s.pop_back(); //saco los :
			int n = stoi(s);
			n++;

			do {
				inputStream >> s;
				if (s == "-1")
					continue;
				outputStream << n << " " << stoi(s)+1 << endl;
			} while (s != "-1");
		}
}

bool nombreOrdenar(pair<int,int> p1,pair<int,int> p2){

	if(p1.first == p2.first)
	{
		return p1.second < p2.second; 
	} else {
		return p1.first < p2.first;
	}
}

void sortEdgesOfGraph(string input_file, string output_file){

		ifstream inputStream(input_file);
		ofstream outputStream(output_file);
		string s;

		int nodes, edges;

		inputStream >> nodes >> edges;

		int from, to;
		vector<pair<int,int> > alledges;
//leer
		while (inputStream >> from) {
			inputStream >> to;
			alledges.push_back(make_pair(from,to));

		}

		sort(alledges.begin(),alledges.end(),nombreOrdenar);

//escribir
		outputStream << nodes << " " << edges << endl;

		for(int i=0;i<alledges.size();i++) {
			outputStream<<alledges[i].first << " "<< alledges[i].second << endl;
		}

}



int main(int argc, char* argv[])
{	
	string input_file = argv[1];
	string output_file = argv[2];

	//modifyInputGraph(input_file, output_file);

	sortEdgesOfGraph(input_file, output_file);



/*

	UniformDist udist(0,1);

	if (argc < 3) {

		cout << "Invalid arguments. use ./generateGraph <nodes> <edges> <outfile>" << endl;
		return 0;

	} else {

		int nodes = stoi(argv[1]);
		int edges = stoi(argv[2]);
		string out_file = argv[3];

		UniformDist udist(1,nodes);

		vector<pair<int,int>> alledges;

		while (alledges.size() < edges) {
			int from = udist.Generate();
			int to = udist.Generate();

			pair<int,int> edge = make_pair(from,to);

			if (from != to && !isPair(edge,alledges))
				alledges.push_back(edge);
		}

		ofstream outputFile(out_file);
		outputFile << nodes << " " << edges << endl; 

		for (int i = 0; i < alledges.size(); i++)
			outputFile << alledges[i].first << " " << alledges[i].second << endl; 
	}

	return 0;
*/

}

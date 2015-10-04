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

int main(int argc, char* argv[])
{	
	///generar muchos c uniformemente..


/*

	UniformDist udist(0,1);

	vector<double> cs;
	for(int i=0;i<100;i++){
		double c =udist.Generate();
		cs.push_back(c);
		cout<<c<<endl;
	}

*/


if(argc<3){

	cout<<"Invalid arguments. use ./generateGraph <nodes> <edges> <outfile>"<<endl;
	return 0;

}else{

	int nodes= stoi(argv[1]);
	int edges= stoi(argv[2]);
	string out_file = argv[3];

	UniformDist udist(1,nodes);

	vector<pair<int,int> > alledges;

	while(alledges.size()<edges){

		int from = udist.Generate();
		int to = udist.Generate();

		pair<int,int> edge = make_pair(from,to);

		if(from != to && !isPair(edge,alledges))
			alledges.push_back(edge);

	}


	ofstream outputFile(out_file);
	outputFile << nodes << " " << edges << endl; 

	for (int i = 0; i < alledges.size(); i++)
		outputFile << alledges[i].first << " " << alledges[i].second << endl; 
	}



	return 0;
}
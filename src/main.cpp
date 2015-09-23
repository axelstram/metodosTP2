
#include <string>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#define PAGE_RANK_METHOD 0
#define ALT_METHOD 1

#define WEB_RANK 0
#define SPORT_RANK 1

using namespace std;


///from test.in
int method;
double c;
int instance_type;
string test_in_file;
string graph_file;
double tolerance;

///from test.txt
int nodes;
int edges;


void load_test_in(){

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
		loadWebGraph();
	}else{
		if(instance_type == SPORT_RANK){
			loadSportGraph();
		}
	}


}





int main(int argc, char* argv[])
{


	test_in_file = argv[1];
	load_test_in();



	cout<< method<<" "<<c<<" "<<instance_type<<" "<<test_in_file<<" "<<graph_file<<" "<<tolerance <<endl;
}

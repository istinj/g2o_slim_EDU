#include <stdio.h>
#include <iostream>
#include <set>
#include <Eigen/Core>
//#include <Eigen/Geometry>

#include "Graph.h"
#include "SparseSolver.h"

using namespace std;

int main(int argc, char const *argv[])
{

	string dataset_path;
	std::vector<string> args(argc);
	if(argc < 2){
		cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << endl;
		return 0;
	}

	for(int i = 1; i < argc; ++i){
		args[i] = argv[i];
		if(args[i] == "-h"){
			cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << endl;
		}
		if(args[i] == "-l"){
			dataset_path = argv[i+1];
			cerr << BOLDWHITE << "Loading world from " << dataset_path << RESET << endl;
			i++;
		}
		else{
			cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << endl;
		}
	}

	optimizer::Graph* graph = new optimizer::Graph(dataset_path);
	optimizer::SparseSolver* solver = new optimizer::SparseSolver(graph->verticesSE3(),
			graph->verticesXYZ(), graph->edgesPosePose(), graph->edgesPosePoint(),
			0, 1000.0);

	for (int i = 0; i < 10; ++i) {
		solver->oneStep();
	}

	solver->updateGraph((*graph));


	delete solver;
	delete graph;

	return 0;
}

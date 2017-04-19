#include <stdio.h>
#include <iostream>
#include <set>
#include <chrono>
#include <Eigen/Core>

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


	std::chrono::high_resolution_clock::time_point t_0, t_1;
	for (int i = 0; i < 10; ++i) {
		t_0 = std::chrono::high_resolution_clock::now();
		solver->oneStep();
		t_1 = std::chrono::high_resolution_clock::now();
		double execution_time = (std::chrono::duration_cast<std::chrono::microseconds>(t_1 - t_0).count() / 1e06);
		std::cerr << BOLDWHITE << "Execution time:\t" << BOLDGREEN << execution_time << std::endl << RESET;
	}

	solver->updateGraph((*graph));


	delete solver;
	delete graph;

	return 0;
}

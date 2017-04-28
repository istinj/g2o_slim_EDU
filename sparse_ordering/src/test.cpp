#include "defs.h"
#include "Graph.h"
#include "SparseOptimizer.h"

int main(int argc, char const *argv[]){
	std::string dataset_path;
	std::vector<std::string> args(argc);
	if(argc < 2){
		std::cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << std::endl;
		return 0;
	}

	for(int i = 1; i < argc; ++i){
		args[i] = argv[i];
		if(args[i] == "-h"){
			std::cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << std::endl;
		}
		if(args[i] == "-l"){
			dataset_path = argv[i+1];
			std::cerr << BOLDWHITE << "Loading world from " << dataset_path << RESET << std::endl;
			i++;
		}
		else{
			std::cerr << YELLOW << "Type -l <path-to-world.g20> to load the world from file" << RESET << std::endl;
		}
	}

	sparse::Graph* graph = new sparse::Graph();
	graph->loadFromG2OFile(dataset_path);

	sparse::SparseOptimizer* optimizer = new sparse::SparseOptimizer();
	optimizer->init(graph->vertices(), graph->edges());
	std::cerr << "Init done" << std::endl;

	std::chrono::high_resolution_clock::time_point t_0, t_1;
	for (int i = 0; i < 10; ++i) {
		t_0 = std::chrono::high_resolution_clock::now();
		optimizer->oneStep();
		t_1 = std::chrono::high_resolution_clock::now();
		double execution_time = (std::chrono::duration_cast<std::chrono::microseconds>(t_1 - t_0).count() / 1e06);
		std::cerr << BOLDWHITE << "Execution time:\t" << BOLDGREEN << execution_time << std::endl << RESET;
	}

	delete optimizer;
	delete graph;
	return 0;
}

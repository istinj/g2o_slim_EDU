#include "defs.h"
#include "Graph.h"
#include "SparseOptimizer.h"

int main(int argc, char const *argv[]){
	sparse::Graph* graph = new sparse::Graph();
//	graph->loadFromG2OFile("../data/world_odom_500p_OL.g2o");
	graph->loadFromG2OFile("../data/world_odom_10kp_OL.g2o");

	sparse::SparseOptimizer* optimizer = new sparse::SparseOptimizer();
	optimizer->init(graph->vertices(), graph->edges());

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

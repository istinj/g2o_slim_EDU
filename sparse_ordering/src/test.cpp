#include <iostream>
#include <Eigen/Core>

#include "defs.h"
#include "Graph.h"
#include "SparseOptimizer.h"

int main(int argc, char const *argv[]){
	sparse::Graph* graph = new sparse::Graph();
	graph->loadFromG2OFile("../data/world_odom_500p_OL.g2o");

	sparse::SparseOptimizer* optimizer = new sparse::SparseOptimizer();
	optimizer->init(graph->vertices(), graph->edges());
	optimizer->oneStep();

	delete optimizer;
	delete graph;
	return 0;
}

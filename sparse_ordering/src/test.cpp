#include <iostream>
#include <Eigen/Core>

#include "defs.h"
#include "Graph.h"

int main(int argc, char const *argv[]){
	sparse::Graph* graph = new sparse::Graph();
	graph->loadFromG2OFile("../data/world_odom_500p_OL.g2o");
	return 0;
}

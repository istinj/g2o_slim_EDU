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
//
//
//  //			map[indices] = new SparseMatrixBlock();
//  //			map[indices]->setIdentity();

//  typedef std::map<std::string, SparseMatrixBlock*> StringMap;
//
//  Matrix6*** memory_map = new Matrix6**[50];
//  for (size_t i = 0; i < 50; ++i){
//    memory_map[i] = new Matrix6*[100];
//  }
//
//  sparse::WorkspaceMap pointers_map;
//  StringMap map3;
//  for (int i = 0; i < 50; ++i){
//    for (int j = 0; j < 100; ++j){
//
//
//      Association indices(i,j);
//      std::cerr << GREEN << indices.first << " " << CYAN <<
//          indices.second << std::endl << RESET;
//
//      Matrix6* block = new Matrix6();
//      block->setIdentity();
//
//      memory_map[i][j] = block;
//
////      map3.insert(std::make_pair(std::to_string(i)+":"+std::to_string(j),block));
//
//      pointers_map.insert(std::make_pair(indices, memory_map[i][j]));
//      *pointers_map.at(indices)->setIdentity();
////      Matrix6 block2 = Matrix6::Identity();
////
////      map2.insert(std::make_pair(indices, block2));
//
////      std::cerr << *map3.at(std::to_string(i)+":"+std::to_string(j)) << std::endl;
//      std::cerr << pointers_map.at(Association(i,j)) << std::endl;
//    }
//  }
//  return 0;
}

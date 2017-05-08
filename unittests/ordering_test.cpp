#include "defs.h"
#include "graph.h"
#include "sparse_optimizer.h"

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
  optimizer->init(graph->vertices(), graph->edges(), PLAIN);

  optimizer->converge();

  optimizer->updateGraph(*graph);

  std::string output_file("../data/world_updated.g2o");
  std::cerr << BOLDWHITE << "\n\nExporting optimized world in " << output_file << RESET << std::endl;
  graph->exportToG2OFile(output_file);

  delete optimizer;
  delete graph;
  return 0;
}

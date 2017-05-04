/*
 * Graph.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include "defs.h"
#include "Vertex.h"
#include "Edge.h"

namespace sparse {

class Graph {
public:
  Graph();
  virtual ~Graph();

  void addVertex(const Vertex& v_);
  void addEdge(const Edge& e_);
  //! TODO: REMOVE Vertex/Edge

  void loadFromG2OFile(const std::string& filename_);

  inline const std::vector<Vertex>& vertices(void) const {return _vertices;};
  inline const std::vector<Edge>& edges(void) const {return _edges;};

protected:
  std::vector<Vertex> _vertices;
  std::vector<Edge> _edges;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* GRAPH_H_ */

/*
 * Vertex.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <Vertex.h>

namespace sparse {

Vertex::Vertex() {
  _id = -1;
  _index = -1;
  _data = Pose::Identity();
  _fixed = false;
}

Vertex::Vertex(const int& id_,
    const int& index_, const Pose& data_){
  _id = id_;
  _index = index_;
  _data = data_;
  _fixed = false;
}

Vertex::~Vertex() {
  // TODO Auto-generated destructor stub
}

void Vertex::setVertex(const int& id_,
    const int& index_, const Pose& data_){
  _id = id_;
  _index = index_;
  _data = data_;
}

void Vertex::setData(const Pose& data_){
  _data = data_;
}

void Vertex::setFixed(const bool fixed_){
  _fixed = fixed_;
}

} /* namespace sparse */

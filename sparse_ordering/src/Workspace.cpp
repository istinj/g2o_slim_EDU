/*
 * Workspace.cpp
 *
 *  Created on: 28/apr/2017
 *      Author: istin
 */

#include <Workspace.h>

namespace sparse {

Workspace::Workspace() {
  dimension = 0;

}

Workspace::~Workspace() {
  dimension = 0;
}

void Workspace::reset(void) {
  if (dimension == 0)
    return;
  for (WorkspaceMap::iterator it = memory_map.begin(); it != memory_map.end(); ++it){
    delete it->second;
  }
}

void Workspace::clear(void){
  if (dimension == 0)
    return;
  for (WorkspaceMap::iterator it = memory_map.begin(); it != memory_map.end(); ++it){
    it->second->setZero();
  }
}

void Workspace::allocate(std::vector<Factor>& factors_){
  for (int i = 0; i < factors_.size(); ++i) {
    int a = std::min(factors_[i].from, factors_[i].to);
    int b = std::max(factors_[i].from, factors_[i].to);

    if(!memory_map.count(Association(a,a))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      memory_map.insert(std::make_pair(Association(a,a), block));
    } else {
      memory_map[Association(a,a)]->setIdentity();
    }

    if(!memory_map.count(Association(b,a))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      memory_map.insert(std::make_pair(Association(b,a), block));
    } else {
      memory_map[Association(a,a)]->setIdentity();
    }

    if(!memory_map.count(Association(a,b))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      memory_map.insert(std::make_pair(Association(a,b), block));
    } else {
      memory_map[Association(a,a)]->setIdentity();
    }

    if(!memory_map.count(Association(b,b))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      memory_map.insert(std::make_pair(Association(b,b), block));
    } else {
      memory_map[Association(a,a)]->setIdentity();
    }
  }

  dimension = memory_map.size();
}


void Workspace::allocate(std::vector<Association>& indices_){
  for (int i = 0; i < indices_.size(); ++i) {
    int& r = indices_[i].first;
    int& c = indices_[i].second;

    if(!memory_map.count(Association(r,c))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      memory_map.insert(std::make_pair(Association(r,c), block));
    } else {
      memory_map[Association(r,c)]->setIdentity();
    }
  }
  dimension = memory_map.size();
}


SparseMatrixBlock& Workspace::operator ()(const Association& indices_) const {
  return *memory_map.at(indices_);
}

SparseMatrixBlock& Workspace::operator ()(const int r_, const int c_) const {
  return *memory_map.at(Association(r_,c_));
}



} /* namespace sparse */

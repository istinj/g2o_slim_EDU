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
  for (WorkspaceMap::iterator it = map.begin(); it != map.end(); ++it){
    delete it->second;
  }
}

void Workspace::clear(void){
  if (dimension == 0)
    return;
  for (WorkspaceMap::iterator it = map.begin(); it != map.end(); ++it){
    it->second->setZero();
  }
}

void Workspace::allocate(std::vector<Factor>& factors_){
  for (int i = 0; i < factors_.size(); ++i) {
    int a = std::min(factors_[i].from, factors_[i].to);
    int b = std::max(factors_[i].from, factors_[i].to);

    if(!map.count(Association(a,a))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      map.insert(std::make_pair(Association(a,a), block));
    } else {
      map[Association(a,a)]->setIdentity();
    }

    if(!map.count(Association(b,a))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      map.insert(std::make_pair(Association(b,a), block));
    } else {
      map[Association(a,a)]->setIdentity();
    }

    if(!map.count(Association(a,b))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      map.insert(std::make_pair(Association(a,b), block));
    } else {
      map[Association(a,a)]->setIdentity();
    }

    if(!map.count(Association(b,b))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      map.insert(std::make_pair(Association(b,b), block));
    } else {
      map[Association(a,a)]->setIdentity();
    }
  }

  dimension = map.size();
}


void Workspace::allocate(std::vector<Association>& indices_){
  for (int i = 0; i < indices_.size(); ++i) {
    int& r = indices_[i].first;
    int& c = indices_[i].second;

    if(!map.count(Association(r,c))){
      SparseMatrixBlock* block = new SparseMatrixBlock();
      block->setIdentity();
      map.insert(std::make_pair(Association(r,c), block));
    } else {
      map[Association(r,c)]->setIdentity();
    }
  }
  dimension = map.size();
}

} /* namespace sparse */

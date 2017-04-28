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

void Workspace::setZero(void){
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
			map[Association(a,a)] = new SparseMatrixBlock();
			map[Association(a,a)]->setIdentity();
		}

		if(!map.count(Association(b,a))){
			map[Association(b,a)] = new SparseMatrixBlock();
			map[Association(b,a)]->setIdentity();
		}

//		if(!_workspace_map.count(Association(a,b))){
//			_workspace_map[Association(a,b)] = new SparseMatrixBlock();
//			_workspace_map[Association(a,b)]->setIdentity();
//		}

		if(!map.count(Association(b,b))){
			map[Association(b,b)] = new SparseMatrixBlock();
			map[Association(b,b)]->setIdentity();
		}
	}
	dimension = map.size();
}

} /* namespace sparse */

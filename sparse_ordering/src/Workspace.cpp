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

		bool verbose = false;

		if (verbose) {
			std::cerr << BOLDGREEN << "Factor(" << factors_[i].from << ", " << factors_[i].to << ")" << std::endl << RESET;
		}

		if(!map.count(Association(a,a))){
			SparseMatrixBlock* block = new SparseMatrixBlock();
			block->setIdentity();
			map.insert(std::make_pair(Association(a,a), block));

//			map[Association(a,a)] = new SparseMatrixBlock();
//			map[Association(a,a)]->setIdentity();
			if (verbose) {
				std::cerr << "Inserted (" << a << ", " << a << ")"<< std::endl;
			}
		} else {
			if (verbose)
				std::cerr << "Skipping (" << a << ", " << a << ")"<< std::endl;
			map[Association(a,a)]->setIdentity();
		}

		if(!map.count(Association(b,a))){
			SparseMatrixBlock* block = new SparseMatrixBlock();
			block->setIdentity();
			map.insert(std::make_pair(Association(b,a), block));

//			map[Association(b,a)] = new SparseMatrixBlock();
//			map[Association(b,a)]->setIdentity();
			if (verbose) {
				std::cerr << "Inserted (" << b << ", " << a << ")"<< std::endl;
			}
		} else {
			if (verbose)
				std::cerr << "Skipping (" << b << ", " << a << ")"<< std::endl;
			map[Association(a,a)]->setIdentity();
		}

/*
		if(!map.count(Association(a,b))){
			map[Association(a,b)] = new SparseMatrixBlock();
			map[Association(a,b)]->setIdentity();
			if (verbose) {
				std::cerr << "Inserted (" << a << ", " << b << ")"<< std::endl;
			}
		} else {
			if (verbose)
				std::cerr << "Skipping (" << a << ", " << b << ")"<< std::endl;
			map[Association(a,a)]->setIdentity();
		}
/**/

		if(!map.count(Association(b,b))){
			SparseMatrixBlock* block = new SparseMatrixBlock();
			block->setIdentity();
			map.insert(std::make_pair(Association(b,b), block));

//			map[Association(b,b)] = new SparseMatrixBlock();
//			map[Association(b,b)]->setIdentity();
			if (verbose) {
				std::cerr << "Inserted (" << b << ", " << b << ")"<< std::endl;
			}
		} else {
			if (verbose)
				std::cerr << "Skipping (" << b << ", " << b << ")"<< std::endl;
			map[Association(a,a)]->setIdentity();
		}
		if (verbose) {
			std::cin.get();
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
			//! WHY THIS SEGFAULTS???
//			map[Association(r,c)] = new SparseMatrixBlock();
//			map[Association(r,c)]->setIdentity();
		} else {
			map[Association(r,c)]->setIdentity();
		}
	}
	dimension = map.size();
}

} /* namespace sparse */

/*
 * SparseBlockMatrix.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <SparseBlockMatrix.h>

using namespace std;

namespace sparse {

SparseBlockMatrix::SparseBlockMatrix() {
	_num_block_rows = 0;
	_num_block_cols = 0;
}

SparseBlockMatrix::SparseBlockMatrix(const VerticesContainer& vertices_,
		const FactorsContainer& factors_){
	_num_block_rows = vertices_.size();
	_num_block_cols = vertices_.size();
	_block_rows.resize(_num_block_rows);

	for (int r = 0; r < _num_block_rows; ++r) {
		for (int c = 0; c < _num_block_cols; ++c) {
			//! TODO: is this a real shit or just a proper shit?
			FactorsContainer::const_iterator factor = find(factors_.begin(),
					factors_.end(), Association(r,c));

		}
	}
}

SparseBlockMatrix::SparseBlockMatrix(const std::vector<Vertex>& vertices_,
		const std::vector<Factor>& factors_,
		const std::vector<int> ordering_){

}

SparseBlockMatrix::~SparseBlockMatrix() {
	// TODO Auto-generated destructor stub
}

} /* namespace sparse */

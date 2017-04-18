/*
 * Factor.h
 *
 *  Created on: 17/apr/2017
 *      Author: istin
 */

#ifndef FACTOR_H_
#define FACTOR_H_

#include "defs.h"

namespace sparse {

struct FactorComparator {
	inline bool operator () (const Factor& a, const Factor& b){
		if(a.from < b.from){
			return true;
		} else {
			if(a.to < b.to)
				return true;
		}
		return false;
	}
};


struct Factor {
	int from;
	int to;

	SparseMatrixBlock* block_ii;
	SparseMatrixBlock* block_ji;
	SparseMatrixBlock* block_jj;

	bool operator==(const Factor& other_) const;
	bool operator==(const Association& factor_association_) const;

	void clear(void);

	//! TODO: use shared ptr
//	std::shared_ptr<SparseMatrixBlock> block_11;
//	std::shared_ptr<SparseMatrixBlock> block_21;
//	std::shared_ptr<SparseMatrixBlock> block_22;
};

} /* namespace sparse */

#endif /* FACTOR_H_ */

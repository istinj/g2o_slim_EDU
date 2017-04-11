/*
 * Factor.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef FACTOR_H_
#define FACTOR_H_

#include "defs.h"

namespace sparse {

class Factor {
public:
	Factor();
	virtual ~Factor();

protected:
	int _from;
	int _to;

	std::vector<SparseMatrixBlock*, Eigen::aligned_allocator<SparseMatrixBlock*> > _factors;

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* FACTOR_H_ */

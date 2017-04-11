/*
 * SparseOptimizer.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef SPARSEOPTIMIZER_H_
#define SPARSEOPTIMIZER_H_

#include "defs.h"

namespace sparse {

class SparseOptimizer {
public:
	SparseOptimizer();
	virtual ~SparseOptimizer();

protected:
	//! TODO: container + SparseBlockMatrix + DenseVector + LinearSolver
};

} /* namespace sparse */

#endif /* SPARSEOPTIMIZER_H_ */

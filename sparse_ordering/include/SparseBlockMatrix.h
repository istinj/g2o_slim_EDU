/*
 * SparseBlockMatrix.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef SPARSEBLOCKMATRIX_H_
#define SPARSEBLOCKMATRIX_H_

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <Eigen/Core>

#include "defs.h"

namespace sparse {
typedef std::map<int, SparseMatrixBlock*, std::less<int>,
		Eigen::aligned_allocator<std::pair<const int, SparseMatrixBlock*> > > ColumnsMap;
typedef std::vector<ColumnsMap> RowsContainer;


class SparseBlockMatrix {
public:
	SparseBlockMatrix();
	virtual ~SparseBlockMatrix();

protected:
	int _total_rows;
	int _total_cols;
	int _num_block_rows;
	int _num_block_cols;

	RowsContainer _block_rows;

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* SPARSEBLOCKMATRIX_H_ */

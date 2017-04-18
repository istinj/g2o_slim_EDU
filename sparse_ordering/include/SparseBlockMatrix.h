/*
 * SparseBlockMatrix.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef SPARSEBLOCKMATRIX_H_
#define SPARSEBLOCKMATRIX_H_

#include "defs.h"
#include "Factor.h"
#include "Vertex.h"

namespace sparse {

typedef std::vector<sparse::Vertex> VerticesContainer;
typedef std::vector<sparse::Factor> FactorsContainer;

typedef std::map<int, SparseMatrixBlock*, std::less<int>,
		Eigen::aligned_allocator<std::pair<const int, SparseMatrixBlock*> > > ColumnsMap;
typedef std::vector<ColumnsMap> RowsContainer;


class SparseBlockMatrix {
public:
	SparseBlockMatrix();
	SparseBlockMatrix(const VerticesContainer& vertices_,
			const FactorsContainer& factors_);
	SparseBlockMatrix(const VerticesContainer& vertices_,
			const FactorsContainer& factors_,
			const std::vector<int> ordering_);
	virtual ~SparseBlockMatrix();

protected:
	int _num_block_rows;
	int _num_block_cols;

	RowsContainer _block_rows;

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* SPARSEBLOCKMATRIX_H_ */

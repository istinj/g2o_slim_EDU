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
typedef std::map<Factor, SparseMatrixBlock*, FactorComparator> FactorsMap;

typedef std::map<int, SparseMatrixBlock*, std::less<int>,
		Eigen::aligned_allocator<std::pair<const int, SparseMatrixBlock*> > > ColumnsMap;
typedef std::vector<ColumnsMap> RowsContainer;


class SparseBlockMatrix {
public:
	SparseBlockMatrix();
	SparseBlockMatrix(const int num_block_rows_,
			const int num_block_cols_, bool has_storage_ = false);
	SparseBlockMatrix(const VerticesContainer& vertices_,
			const FactorsMap& factors_);
	SparseBlockMatrix(const VerticesContainer& vertices_,
			const FactorsMap& factors_,
			const std::vector<int> ordering_);
	virtual ~SparseBlockMatrix();

	inline const int numRows(void) const {return _num_block_rows;};
	inline const int numCols(void) const {return _num_block_cols;};

	void setBlock(const int r_,
			const int c_,
			SparseMatrixBlock* data_ptr_);
	SparseMatrixBlock getBlock(const int r_, const int c_) const;
	SparseMatrixBlock* getBlockPtr(const int r_, const int c_) const;

	SparseBlockMatrix* transpose(void) const;
	//!TODO: 	this is a shit but does not produce memory leaks
	SparseBlockMatrix* rightMultiplySparseMatrix(const SparseBlockMatrix* other_) const;
	//!TODO: 	this is a shit but does not produce memory leaks
	SparseBlockMatrix* cholesky(void) const;

	void printBlock(const int r_, const int c_) const;
	void printMatrix(void) const;
protected:
	SparseMatrixBlock* scalarProdPtr(const ColumnsMap& row1_,
			const ColumnsMap& row2_,
			const int max_pos_) const;
	//!TODO: 	this is a shit but does not produce memory leaks
	std::shared_ptr<SparseMatrixBlock> scalarProd(const ColumnsMap& row1_,
			const ColumnsMap& row2_,
			const int max_pos_) const;

	int _num_block_rows;
	int _num_block_cols;
	bool _has_storage;

	RowsContainer _block_rows;
	std::map<Association, SparseMatrixBlock*, AssociationComparator> _storage;

	//! TODO 	MEMORY MANAGEMENT when the matrix owns the blocks.
	//! TODO 	ORDERING

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* SPARSEBLOCKMATRIX_H_ */

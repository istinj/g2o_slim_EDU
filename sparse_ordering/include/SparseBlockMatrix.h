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

struct DenseBlockVector {
	int size;
	DenseVectorContainer blocks;

	inline void reset(void){
		for (int i = 0; i < size; ++i) {
			delete blocks[i];
		}
	}

	inline void clear(void) {
		for (int i = 0; i < blocks.size(); ++i) {
			blocks[i]->setZero();
		}
	}

	inline void init(const int size_){
		size = size_;
		for (int i = 0; i < size; ++i) {
			blocks[i] = new DenseVectorBlock();
			blocks[i]->setZero();
		}
	}
};

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

	void printBlock(const int r_, const int c_) const;
	void printMatrix(void) const;

	void setBlock(const int r_,	const int c_, SparseMatrixBlock* data_ptr_);
	SparseMatrixBlock getBlock(const int r_, const int c_) const;
	SparseMatrixBlock* getBlockPtr(const int r_, const int c_) const;

	SparseBlockMatrix* transpose(void) const;
	//! TODO:	This shit segfaulta
	SparseBlockMatrix* rightMultiplySparseMatrix(const SparseBlockMatrix* other_) const;
	SparseBlockMatrix* cholesky(void) const;

	//! This produces memory access.
	void solveLinearSystem(DenseBlockVector& rhs_vector_, DenseBlockVector& result_) const;
	void forwSub(DenseBlockVector& rhs_vector_, DenseBlockVector& result_) const;
	void backSub(DenseBlockVector& rhs_vector_, DenseBlockVector& result_) const;

protected:
	//!TODO: 	this is a shit but does not produce memory leaks
	SparseMatrixBlock scalarProd(const ColumnsMap& row1_,
			const ColumnsMap& row2_,
			const int max_pos_) const;

	int _num_block_rows;
	int _num_block_cols;
	bool _has_storage;

	RowsContainer _block_rows;
	std::map<Association, SparseMatrixBlock*, AssociationComparator> _storage;


	//! TODO	Method to remove row/column
	//! TODO	Operator() overload
	//! TODO	CHECK MEMORY MANAGEMENT
	//! TODO 	MEMORY MANAGEMENT when the matrix owns the blocks.
	//! TODO 	ORDERING

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* SPARSEBLOCKMATRIX_H_ */

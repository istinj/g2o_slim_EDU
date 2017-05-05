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
#include "Workspace.h"

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

typedef std::vector<sparse::Vertex> VerticesVector;

typedef std::map<int, SparseMatrixBlock*, std::less<int>,
        Eigen::aligned_allocator<std::pair<const int, SparseMatrixBlock*> > > IntSparseMatrixBlockPtrMap;
typedef std::vector<IntSparseMatrixBlockPtrMap> IntSparseMatrixBlockPtrMapVector;
typedef std::vector<Factor> FactorsVector;


class SparseBlockMatrix {
public:
  SparseBlockMatrix();
  SparseBlockMatrix(const int num_block_rows_,
                    const int num_block_cols_,
                    bool has_storage_ = false);
  SparseBlockMatrix(const int size_,
                    const Workspace& workspace_);
  virtual ~SparseBlockMatrix();

  inline const int numRows(void) const {return _num_block_rows;}
  inline const int numCols(void) const {return _num_block_cols;}
  inline const int nnz(void) const {return _nnz;}

  void setZero(void);
  void printBlock(const int r_, const int c_) const;
  void printMatrix(void) const;

  void setBlock(const int r_,
                const int c_,
                SparseMatrixBlock* data_ptr_);
  SparseMatrixBlock getBlock(const int r_,
                             const int c_) const;
  SparseMatrixBlock* getBlockPtr(const int r_,
                                 const int c_) const;
  bool isNonZeroBlock(const int r_,
                      const int c_) const;

  void reorder(const IntVector& permutation_,
               const Workspace& workspace_,
               const FactorsVector& factors_);

  void allocateTransposed(SparseBlockMatrix& transposed_);
  void computeTranspose(SparseBlockMatrix& transposed_);

  void allocateCholesky(SparseBlockMatrix& cholesky_);
  void computeCholesky(SparseBlockMatrix& cholesky_);

  void forwSub(DenseBlockVector& rhs_vector_,
               DenseBlockVector& result_) const;
  void backSub(DenseBlockVector& rhs_vector_,
               DenseBlockVector& result_) const;

  cs* toCs(void) const;

protected:
  SparseMatrixBlock scalarProd(const IntSparseMatrixBlockPtrMap& row1_,
                               const IntSparseMatrixBlockPtrMap& row2_,
                               const int max_pos_) const;
  bool scalarProdStructure(const IntSparseMatrixBlockPtrMap& row1_,
                           const IntSparseMatrixBlockPtrMap& row2_,
                           const int max_pos_) const;

  int _num_block_rows;
  int _num_block_cols;
  int _nnz;
  bool _has_storage;
  bool _is_initialized = false;

  IntSparseMatrixBlockPtrMapVector _block_rows;
  Workspace _matrix_workspace;


  //! TODO 	ORDERING
  //! TODO	Method to remove row/column
  //! TODO	Operator() overload

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

std::ostream& operator <<(std::ostream& os, const SparseBlockMatrix& m);

} /* namespace sparse */

#endif /* SPARSEBLOCKMATRIX_H_ */

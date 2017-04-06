/*
 * SparseBlockMatrix.hpp
 *
 *  Created on: 29/mar/2017
 *      Author: istin
 */

#ifndef SPARSEBLOCKMATRIX_HPP_
#define SPARSEBLOCKMATRIX_HPP_

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <boost/unordered_map.hpp>

#include "utilities.h"
#include "RHSVector.h"

namespace sparse {

template<typename BlockType_>
class SparseBlockMatrix {
public:
   typedef BlockType_ DenseBlock; //! TODO Must be a static Eigen::Matrix -> manage exceptions
   typedef std::map<int, DenseBlock, std::less<int>,
         Eigen::aligned_allocator<std::pair<const int, DenseBlock> > > ColumnsBlockMap;

   typedef std::vector<ColumnsBlockMap> RowBlockContainer;

   SparseBlockMatrix();
   SparseBlockMatrix(const int num_rows_, const int num_cols_, const int block_dim_);
   virtual ~SparseBlockMatrix();

   void resize(const int new_rows_, const int new_cols_);
   void reset(void);

   void computeCholeskyStructure(SparseBlockMatrix<BlockType_>& block_cholesky_);
   bool computeScalarProdStructure(const ColumnsBlockMap& row_1_,
         const ColumnsBlockMap& row_2_, int max_pos_);

   void transpose(SparseBlockMatrix<BlockType_>& transpose_);
   void cholesky(SparseBlockMatrix<BlockType_>& block_cholesky_);
   DenseBlock scalarProd(const ColumnsBlockMap& row_1_,
         const ColumnsBlockMap& row_2_, int max_pos_);
   //! TODO: pass a pointer?
   template<typename VectorBlockType_>
   bool solveLinearSystem(const DenseVector<VectorBlockType_>& RHS_Vector_,
         DenseVector<VectorBlockType_>& result_);
   template<typename VectorBlockType_>
   void forwSubstitution(const DenseVector<VectorBlockType_>& B_vector_,
         DenseVector<VectorBlockType_>& result_);
   template<typename VectorBlockType_>
   void backSubstitution(const DenseVector<VectorBlockType_>& B_vector_,
         DenseVector<VectorBlockType_>& result_);


   void setBlock(const int r_, const int c_, DenseBlock block_);
   DenseBlock getBlock(const int r_, const int c_) const;
   bool isNonZeroBlock(const int r_, const int c_) const;
   void printBlock(const int r_, const int c_) const;
   void printMatrix(void) const;

   inline const int numRows(void) const {return _num_block_rows;};
   inline const int numCols(void) const {return _num_block_cols;};

protected:
   int _total_rows;
   int _total_cols;
   int _num_block_rows;
   int _num_block_cols;
   int _block_dim = 1;
   //! TODO: traits to manage wrong BlockType_
   //! TODO: handle non-squared blocks (pair in the constructor);

   RowBlockContainer _row_container;

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#include "SparseBlockMatrix.hpp"
#endif /* SPARSEBLOCKMATRIX_HPP_ */

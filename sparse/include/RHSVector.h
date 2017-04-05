/*
 * RHSVector.h
 *
 *  Created on: 03/apr/2017
 *      Author: istin
 */

#ifndef RHSVECTOR_H_
#define RHSVECTOR_H_

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

namespace sparse {

template<typename VectorBlockType_>
class DenseVector {
public:
   typedef VectorBlockType_ DenseVectorBlock; //! TODO Must be a static Eigen::Matrix -> manage exceptions
   typedef std::vector<DenseVectorBlock, Eigen::aligned_allocator<DenseVectorBlock> > DenseVectorContainer;

   DenseVector();
   DenseVector(const int num_block_rows_, const int block_dim_);
   virtual ~DenseVector();

   void resize(const int new_rows_);
   void reset(void);

   void setBlock(const int r_, DenseVectorBlock block_);
   DenseVectorBlock getBlock(const int r_) const;
   void printBlock(const int r_) const;
   void printVector(void) const;

   inline const int numRows(void) const {return _num_block_rows;};

private:
   int _num_block_rows;
   int _total_rows;
   int _block_dim = 1;
   DenseVectorContainer _rhs_block_container;

public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#include "RHSVector.hpp"
#endif /* RHSVECTOR_H_ */

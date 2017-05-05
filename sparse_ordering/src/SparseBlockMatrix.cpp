/*
 * SparseBlockMatrix.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <SparseBlockMatrix.h>

using namespace std;

namespace sparse {

SparseBlockMatrix::SparseBlockMatrix() {
  _num_block_rows = 0;
  _num_block_cols = 0;
  _nnz = 0;
  _has_storage = false;
}

SparseBlockMatrix::SparseBlockMatrix(const int num_block_rows_,
                                     const int num_block_cols_,
                                     bool has_storage_) {
  _num_block_rows = num_block_rows_;
  _num_block_cols = num_block_cols_;
  _block_rows.resize(_num_block_rows);
  _has_storage = has_storage_;
  _is_initialized = true;
  _nnz = 0;
}


SparseBlockMatrix::SparseBlockMatrix(const int size_,
                                     const Workspace& workspace_){
  _num_block_rows = size_;
  _num_block_cols = size_;
  _block_rows.resize(_num_block_rows);

  for (WorkspaceMap::const_iterator it = workspace_.memory_map.begin(); it != workspace_.memory_map.end(); ++it){
    setBlock(it->first.first, it->first.second, it->second);
  }

  _has_storage = false;
  _is_initialized = true;
  _nnz = workspace_.dimension;
}

SparseBlockMatrix::~SparseBlockMatrix() {
  if (_has_storage)
    _matrix_workspace.reset();
  _is_initialized = false;
}

void SparseBlockMatrix::setZero(void) {
  if (_has_storage) {
    _matrix_workspace.clear();
  }
}

void SparseBlockMatrix::setBlock(const int r_,
                                 const int c_,
                                 SparseMatrixBlock* data_ptr_) {
  if(r_ >= _num_block_rows || c_ >= _num_block_cols){
    throw std::runtime_error("You are trying to set an wrong element - out of bound index");
  }

  IntSparseMatrixBlockPtrMap& block_row = _block_rows[r_];
  IntSparseMatrixBlockPtrMap::iterator block = block_row.find(c_);

  if (block == block_row.end()){
    if(!data_ptr_->isZero())
      block_row[c_] = data_ptr_;
  } else {
    if(!data_ptr_->isZero())
      block_row.erase(block);
    else
      block->second = data_ptr_;
  }
}

SparseMatrixBlock SparseBlockMatrix::getBlock(const int r_,
                                              const int c_) const {
  SparseMatrixBlock block = SparseMatrixBlock::Zero();
  if(r_ >= _num_block_rows || c_ >= _num_block_cols){
    throw std::runtime_error("You are trying to access an wrong element - out of bound index");
  }
  const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r_];
  IntSparseMatrixBlockPtrMap::const_iterator it = block_row.find(c_);
  if (it == block_row.end())
    return block;
  return (*it->second);
}

SparseMatrixBlock* SparseBlockMatrix::getBlockPtr(const int r_,
                                                  const int c_) const {
  if(r_ >= _num_block_rows || c_ >= _num_block_cols){
    throw std::runtime_error("You are trying to access an wrong element - out of bound index");
  }

  const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r_];
  IntSparseMatrixBlockPtrMap::const_iterator it = block_row.find(c_);
  if (it == block_row.end())
    return nullptr;
  return it->second;

}

bool SparseBlockMatrix::isNonZeroBlock(const int r_,
                                       const int c_) const {
  const SparseMatrixBlock& block = getBlock(r_, c_);
  if(!block.isZero())
    return true;
  return false;
}


void SparseBlockMatrix::reorder(const IntVector& permutation_,
                                const Workspace& workspace_,
                                const FactorsVector& factors_){
  //! Cleaning up the matrix
  for (int i = 0; i < _num_block_rows; ++i) {
    _block_rows[i].clear();
  }
  _block_rows.clear();
  _block_rows.resize(_num_block_rows);

//  for (size_t i = 0; i < permutation_.size(); ++i) {
//    cerr << permutation_[i] << endl;
//  }

  //! Set again the block-pointers to the right memory blocks
  for (int i = 0; i < factors_.size(); ++i) {
    const Factor& factor = factors_[i];
    int from_idx = permutation_[factor.from];
    int to_idx = permutation_[factor.to];

    int a = std::min(from_idx, to_idx);
    int b = std::max(from_idx, to_idx);

    int old_a = std::min(factor.from, factor.to);
    int old_b = std::max(factor.from, factor.to);

    setBlock(a,a,workspace_.memory_map.at(Association(old_a,old_a)));
    setBlock(a,b,workspace_.memory_map.at(Association(old_a,old_b)));
    setBlock(b,a,workspace_.memory_map.at(Association(old_b,old_a)));
    setBlock(b,b,workspace_.memory_map.at(Association(old_b,old_b)));
  }
}



void SparseBlockMatrix::allocateCholesky(SparseBlockMatrix& cholesky_) {
  if(_num_block_rows != _num_block_cols){
    throw std::runtime_error("Error: matrix must be squared");
  }

  cholesky_ = SparseBlockMatrix(_num_block_rows, _num_block_cols, true);
  std::vector<Association> indices;

  int nnz = 0;
  for (int r = 0; r < _num_block_rows; ++r) {
    const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r];
    IntSparseMatrixBlockPtrMap& chol_block_row = cholesky_._block_rows[r];

    if(block_row.empty())
      throw std::runtime_error("Something went wrong during the cholesky allocation");
    int starting_col_idx = block_row.begin()->first;

    for (int c = starting_col_idx; c <= r; ++c) {
      bool not_empty = false;
      if(isNonZeroBlock(r,c)){
        not_empty = true;
      } else {
        IntSparseMatrixBlockPtrMap& chol_upper_row = cholesky_._block_rows[c];
        not_empty = scalarProdStructure(chol_block_row, chol_upper_row, c-1);
      }
      if(not_empty) {
        indices.push_back(Association(r,c));
        cholesky_._matrix_workspace.allocateOneBlock(r,c);
        cholesky_.setBlock(r,c,cholesky_._matrix_workspace.memory_map[Association(r,c)]);
        ++nnz;
      }
    }
  }
  cholesky_._nnz = nnz;
}


void SparseBlockMatrix::computeCholesky(SparseBlockMatrix& cholesky_) {
  if (!cholesky_._is_initialized)
    throw std::runtime_error("Argument matrix must be already initialized");

  std::vector<SparseMatrixBlock, Eigen::aligned_allocator<SparseMatrixBlock> > inverse_diag_blocks(_num_block_rows);

  for (int r = 0; r < _num_block_rows; ++r) {
    const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r];
    IntSparseMatrixBlockPtrMap& chol_block_row = cholesky_._block_rows[r];

    if(block_row.empty())
      throw std::runtime_error("Something went wrong during the cholesky computation");
    int starting_col_idx = block_row.begin()->first;

    for (int c = starting_col_idx; c <= r; ++c) {

      //! If this block has not been allocated then skip.
      if(cholesky_._matrix_workspace.memory_map[Association(r,c)] == nullptr)
        continue;

      SparseMatrixBlock& chol_computed_block = cholesky_._matrix_workspace(r,c);
      SparseMatrixBlock accumulator = SparseMatrixBlock::Zero();
      chol_computed_block.setZero();

      const IntSparseMatrixBlockPtrMap& chol_upper_row = cholesky_._block_rows[c];
      accumulator = getBlock(r,c) - scalarProd(chol_block_row, chol_upper_row, c-1);
      if (r == c) {
        chol_computed_block = accumulator.llt().matrixL();
        inverse_diag_blocks[r] = chol_computed_block.inverse().transpose();
      } else {
        chol_computed_block = accumulator * inverse_diag_blocks[c];
      }
      chol_block_row[c] = cholesky_._matrix_workspace.memory_map.at(Association(r,c));
    }
  }
}


void SparseBlockMatrix::allocateTransposed(SparseBlockMatrix& transposed_) {
  transposed_ = SparseBlockMatrix(_num_block_rows, _num_block_cols, true);
  std::vector<Association> indices;

  for (int r = 0; r < _num_block_rows; ++r) {
    const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r];
    if(block_row.empty())
      continue;
//      throw std::runtime_error("Something went wrong during the allocation of the transpose");
    for (int c = 0; c < _num_block_cols; ++c) {
      IntSparseMatrixBlockPtrMap& transposed_block_row = transposed_._block_rows[c];

      IntSparseMatrixBlockPtrMap::const_iterator it = block_row.find(c);
      if (it == block_row.end())
        continue;
      else{
        indices.push_back(Association(c,r));
      }
    }
  }
  transposed_._matrix_workspace.allocate(indices);
  for(WorkspaceMap::const_iterator it = transposed_._matrix_workspace.memory_map.begin(); it != transposed_._matrix_workspace.memory_map.end();++it){
    transposed_.setBlock(it->first.first, it->first.second, it->second);
  }
}

void SparseBlockMatrix::computeTranspose(SparseBlockMatrix& transposed_) {
  if (!transposed_._is_initialized)
    throw std::runtime_error("Argument matrix must be already initialized");

  for (int r = 0; r < _num_block_rows; ++r) {
    const IntSparseMatrixBlockPtrMap& block_row = _block_rows[r];
    if(block_row.empty())
      continue;
//      throw std::runtime_error("Something went wrong during the transpose update");
    for (int c = 0; c < _num_block_cols; ++c) {
      IntSparseMatrixBlockPtrMap& transposed_block_row = transposed_._block_rows[c];

      IntSparseMatrixBlockPtrMap::const_iterator it = block_row.find(c);
      if (it == block_row.end())
        continue;
      else{
        transposed_._matrix_workspace(c,r) = it->second->transpose();
        transposed_block_row[r] = transposed_._matrix_workspace.memory_map[Association(c,r)];
      }
    }
  }
}


void SparseBlockMatrix::forwSub(DenseBlockVector& rhs_vector_,
                                DenseBlockVector& result_) const {
  if(_num_block_rows != _num_block_cols || rhs_vector_.size != _num_block_rows)
    throw std::runtime_error("Error, dimensions must agree");

  if (result_.size != rhs_vector_.size) {
    result_.reset();
    result_.init(rhs_vector_.size);
  }

  for (int r = 0; r < _num_block_rows; ++r) {
    DenseVectorBlock& res_block = (*result_.blocks[r]);
    res_block = (*rhs_vector_.blocks[r]);
    for (int c = 0; c < r; ++c) {
      res_block.noalias() -= getBlock(r,c) * (*result_.blocks[c]);
    }
    res_block = getBlock(r,r).inverse() * res_block;
  }
}


void SparseBlockMatrix::backSub(DenseBlockVector& rhs_vector_,
                                DenseBlockVector& result_) const {
  if(_num_block_rows != _num_block_cols || rhs_vector_.size != _num_block_rows)
    throw std::runtime_error("Error, dimensions must agree");

  if (result_.size != rhs_vector_.size) {
    result_.reset();
    result_.init(rhs_vector_.size);
  }

  for (int r = _num_block_rows - 1; r >= 0; --r) {
    DenseVectorBlock& res_block = (*result_.blocks[r]);
    res_block = (*rhs_vector_.blocks[r]);
    for (int c = r + 1; c < _num_block_cols; ++c) {
      res_block.noalias() -= getBlock(r,c) * (*result_.blocks[c]);
    }
    res_block = getBlock(r,r).inverse() * res_block;
  }
}


cs* SparseBlockMatrix::toCs(void) const {
  cs* transposed_matrix = cs_spalloc(_num_block_cols, _num_block_rows, _nnz, 1, 0);

  int *row_pointers = transposed_matrix->p;
  int* col_indices = transposed_matrix->i;
  double* values = transposed_matrix->x;

  int counter=0;
  for (size_t r = 0; r < _num_block_rows; ++r){
    *row_pointers = counter;
    const IntSparseMatrixBlockPtrMap& row = _block_rows[r];
    for (IntSparseMatrixBlockPtrMap::const_iterator it = row.begin(); it != row.end(); ++it) {
      int c = it->first;
      float value = 1.0;
      *values = value;
      *col_indices = c;
      ++counter;
      ++values;
      ++col_indices;
    }
    ++row_pointers;
  }
  *row_pointers = counter;
  cs* tm = cs_transpose(transposed_matrix, 1);
  cs_spfree(transposed_matrix);

  return tm;
}


void SparseBlockMatrix::printBlock(const int r_, const int c_) const {
  if(r_ >= _num_block_rows || c_ >= _num_block_cols){
    throw std::runtime_error("You are trying to print an wrong element - out of bound index");
  }

  cerr << BOLDWHITE << "Block(" << r_ << "," << c_ << ")" << ":\n" <<
      CYAN << getBlock(r_,c_) << RESET << endl;
}

void SparseBlockMatrix::printMatrix(void) const {
  cerr << BOLDWHITE <<  "Printing sparse matrix..." << RESET << endl;
  for (int i = 0; i < _num_block_rows; ++i) {
    for (int j = 0; j < _num_block_cols; ++j)
      if(isNonZeroBlock(i,j))
        printBlock(i,j);
  }
}

SparseMatrixBlock SparseBlockMatrix::scalarProd(const IntSparseMatrixBlockPtrMap& row1_,
                                                const IntSparseMatrixBlockPtrMap& row2_,
                                                const int max_pos_) const {

  typename IntSparseMatrixBlockPtrMap::const_iterator it_1 = row1_.begin();
  typename IntSparseMatrixBlockPtrMap::const_iterator it_2 = row2_.begin();

  SparseMatrixBlock result = SparseMatrixBlock::Zero();
  while(it_1 != row1_.end() && it_2 != row2_.end()){
    int col_idx_1 = it_1->first;
    int col_idx_2 = it_2->first;
    if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_){
      return result;
    }
    if(col_idx_1 == col_idx_2){
      result += (*it_1->second) * it_2->second->transpose();
      ++it_1;
      ++it_2;
    }
    else if(col_idx_1 > col_idx_2)
      ++it_2;
    else if(col_idx_1 < col_idx_2)
      ++it_1;
  }
  return result;
}

bool SparseBlockMatrix::scalarProdStructure(const IntSparseMatrixBlockPtrMap& row1_,
                                            const IntSparseMatrixBlockPtrMap& row2_,
                                            const int max_pos_) const {

  typename IntSparseMatrixBlockPtrMap::const_iterator it_1 = row1_.begin();
  typename IntSparseMatrixBlockPtrMap::const_iterator it_2 = row2_.begin();
  while(it_1 != row1_.end() && it_2 != row2_.end()){
    int col_idx_1 = it_1->first;
    int col_idx_2 = it_2->first;
    if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_)
      return false;
    if(col_idx_1 == col_idx_2)
      return true;
    else if(col_idx_1 > col_idx_2)
      ++it_2;
    else if(col_idx_1 < col_idx_2)
      ++it_1;
  }
  return false;
}


std::ostream& operator <<(std::ostream& os, const SparseBlockMatrix& m){
//  os << m.numRows() << " " << m.numCols() << std::endl;
  for (int r=0; r<m.numRows(); r++){
    for (int c=0; c<m.numCols(); c++){
      if(m.isNonZeroBlock(r,c))
        os << 1 << " ";
      else
        os << 0 << " ";
    }
    os << std::endl;
  }
  return os;
}






} /* namespace sparse */

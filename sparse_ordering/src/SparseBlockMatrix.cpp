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
}

SparseBlockMatrix::SparseBlockMatrix(const VerticesContainer& vertices_,
		const BlocksMap& blocks_){
	_num_block_rows = vertices_.size();
	_num_block_cols = vertices_.size();
	_block_rows.resize(_num_block_rows);

	for (BlocksMap::const_iterator block = blocks_.begin(); block != blocks_.end(); ++block) {
		setBlock(block->first.first, block->first.second, block->second);
	}
	_has_storage = false;
	_is_initialized = true;
}

SparseBlockMatrix::SparseBlockMatrix(const std::vector<Vertex>& vertices_,
		const BlocksMap& blocks_,
		const FactorsVector& factors_,
		const std::vector<int> ordering_){
	_num_block_rows = vertices_.size();
	_num_block_cols = vertices_.size();
	_block_rows.resize(_num_block_rows);

	//! TODO ORDERING (???) - shit over shit
	for (int i = 0; i < factors_.size(); ++i) {
		const Factor& factor = factors_[i];
		int from_idx = factor.from;
		int to_idx = factor.to;

		from_idx = ordering_[from_idx];
		to_idx = ordering_[to_idx];

		int a = std::min(from_idx, to_idx);
		int b = std::max(from_idx, to_idx);

		BlocksMap::const_iterator block_ii = blocks_.find(Association(factor.from, factor.from));
		if(getBlock(a,a) == SparseMatrixBlock::Zero()){
			setBlock(a,a,block_ii->second);
		}

		BlocksMap::const_iterator block_ji = blocks_.find(Association(factor.to, factor.from));
		if(getBlock(b,a) == SparseMatrixBlock::Zero()){
			setBlock(b,a,block_ji->second);
		}

		BlocksMap::const_iterator block_jj = blocks_.find(Association(factor.to, factor.to));
		if(getBlock(b,a) == SparseMatrixBlock::Zero()){
			setBlock(b,b,block_jj->second);
		}

	}

	_has_storage = false;
	_is_initialized = true;
}

SparseBlockMatrix::SparseBlockMatrix(const int size_,
		const Workspace& workspace_){
	_num_block_rows = size_;
	_num_block_cols = size_;
	_block_rows.resize(_num_block_rows);

	//! TODO
	for (WorkspaceMap::const_iterator it = workspace_.map.begin(); it != workspace_.map.end(); ++it){
		setBlock(it->first.first, it->first.second, it->second);
	}

	_has_storage = false;
	_is_initialized = true;
}

//! TODO: how to delete this shit when it owns the memory?
SparseBlockMatrix::~SparseBlockMatrix() {
	if (_has_storage){
		for (std::map<Association, SparseMatrixBlock*, AssociationComparator>::iterator it = _storage.begin(); it != _storage.end(); ++it){
			delete it->second;

			_matrix_workspace.reset();
		}
	}
	_is_initialized = false;
}

//! Why this does not work??
void SparseBlockMatrix::clear(void) {
	if (_has_storage) {
		for (int r = 0; r < 0; ++r) {
			ColumnsMap& block_row = _block_rows[r];
			for(ColumnsMap::iterator it = block_row.begin(); it != block_row.end(); ++it) {
				it->second->setZero();
			}
		}
		_matrix_workspace.setZero();
	}
}

void SparseBlockMatrix::setBlock(const int r_, const int c_, SparseMatrixBlock* data_ptr_) {
	if(r_ >= _num_block_rows || c_ >= _num_block_cols){
		throw std::runtime_error("You are trying to set an wrong element - out of bound index");
	}

	ColumnsMap& block_row = _block_rows[r_];
	ColumnsMap::iterator block = block_row.find(c_);

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

SparseMatrixBlock SparseBlockMatrix::getBlock(const int r_, const int c_) const {
	SparseMatrixBlock block = SparseMatrixBlock::Zero();
	if(r_ >= _num_block_rows || c_ >= _num_block_cols){
		throw std::runtime_error("You are trying to access an wrong element - out of bound index");
	}

	const ColumnsMap& block_row = _block_rows[r_];
	ColumnsMap::const_iterator it = block_row.find(c_);
	if (it == block_row.end())
		return block;
	return (*it->second);

}

SparseMatrixBlock* SparseBlockMatrix::getBlockPtr(const int r_, const int c_) const {
	if(r_ >= _num_block_rows || c_ >= _num_block_cols){
		throw std::runtime_error("You are trying to access an wrong element - out of bound index");
	}

	const ColumnsMap& block_row = _block_rows[r_];
	ColumnsMap::const_iterator it = block_row.find(c_);
	if (it == block_row.end())
		return NULL;
	return it->second;

}

SparseBlockMatrix* SparseBlockMatrix::cholesky(void) const {
	if(_num_block_rows != _num_block_cols){
		throw std::runtime_error("Error: matrix must be squared");
	}
	cerr << RED << _num_block_rows << " " << _num_block_cols << endl << RESET;

	SparseBlockMatrix* chol = new SparseBlockMatrix(_num_block_rows, _num_block_cols, true);
	std::vector<SparseMatrixBlock, Eigen::aligned_allocator<SparseMatrixBlock> > inverse_diag_blocks(_num_block_rows);

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		ColumnsMap& chol_block_row = chol->_block_rows[r];

		if(block_row.empty())
			throw std::runtime_error("Something went wrong during the cholesky");
		int starting_col_idx = block_row.begin()->first;

		for (int c = starting_col_idx; c <= r; ++c) {
			cerr << BOLDGREEN << "(" << r << ", " << c << ")" << endl << RESET;
			cerr << "1" << endl;

			//! Allocate memory
			SparseMatrixBlock* c_block = new SparseMatrixBlock();
			c_block->setZero();
			cerr << "2" << endl;

			chol->_storage.insert(make_pair(Association(r,c), c_block));
			SparseMatrixBlock& chol_computed_block = *c_block;
			cerr << "3" << endl;

//			chol->_storage[Association(r,c)] = new SparseMatrixBlock();
//			cerr << "2" << endl;
//			SparseMatrixBlock& chol_computed_block = (*chol->_storage[Association(r,c)]);
//			chol->_storage[Association(r,c)]->setZero();
////			chol_computed_block.setZero();
//			cerr << "3" << endl;

			SparseMatrixBlock accumulator = SparseMatrixBlock::Zero();

			const ColumnsMap& chol_upper_row = chol->_block_rows[c];
			cerr << BLUE << scalarProd(chol_block_row, chol_upper_row, c-1) << endl << RESET;
			accumulator = getBlock(r,c) - scalarProd(chol_block_row, chol_upper_row, c-1);
			cerr << "4" << endl;

			if (r == c) {
				chol_computed_block = accumulator.llt().matrixL();
				inverse_diag_blocks[r] = chol_computed_block.inverse().transpose();
				cerr << "5a" << endl;
			} else {
				chol_computed_block = accumulator * inverse_diag_blocks[c];
				cerr << "5b" << endl;
			}
			cerr << "6" << endl;
			chol_block_row[c] = chol->_storage[Association(r,c)];
		}
	}
	cerr << "end" << endl;
	return chol;
}

void SparseBlockMatrix::updateCholesky(SparseBlockMatrix* result_) {
	if (!result_->_is_initialized)
		throw std::runtime_error("Argument matrix must be already initialized");

	std::vector<SparseMatrixBlock, Eigen::aligned_allocator<SparseMatrixBlock> > inverse_diag_blocks(_num_block_rows);

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		ColumnsMap& chol_block_row = result_->_block_rows[r];

		if(block_row.empty())
			throw std::runtime_error("Something went wrong during the cholesky update");
		int starting_col_idx = block_row.begin()->first;

		for (int c = starting_col_idx; c <= r; ++c) {
			SparseMatrixBlock& chol_computed_block = (*result_->_storage[Association(r,c)]);

			SparseMatrixBlock accumulator = SparseMatrixBlock::Zero();
			chol_computed_block.setZero();

			const ColumnsMap& chol_upper_row = result_->_block_rows[c];
			accumulator = getBlock(r,c) - scalarProd(chol_block_row, chol_upper_row, c-1);
			if (r == c) {
				chol_computed_block = accumulator.llt().matrixL();
				inverse_diag_blocks[r] = chol_computed_block.inverse().transpose();
			} else {
				chol_computed_block = accumulator * inverse_diag_blocks[c];
			}
			chol_block_row[c] = result_->_storage[Association(r,c)];
		}
	}
}

SparseBlockMatrix* SparseBlockMatrix::transpose(void) const {
	SparseBlockMatrix* transposed = new SparseBlockMatrix(_num_block_rows, _num_block_cols, true);

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		if(block_row.empty())
			throw std::runtime_error("Something went wrong during the transpose");
		for (int c = 0; c < _num_block_cols; ++c) {
			ColumnsMap& transposed_block_row = transposed->_block_rows[c];

			ColumnsMap::const_iterator it = block_row.find(c);
			if (it == block_row.end())
				continue;
			else{
				transposed->_storage[Association(c,r)] = new SparseMatrixBlock();
				(*transposed->_storage[Association(c,r)]) = it->second->transpose();
				transposed_block_row[r] = transposed->_storage[Association(c,r)];
			}
		}
	}
	return transposed;
}

void SparseBlockMatrix::updateTranspose(SparseBlockMatrix* result_) {
	if (!result_->_is_initialized)
		throw std::runtime_error("Argument matrix must be already initialized");

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		if(block_row.empty())
			throw std::runtime_error("Something went wrong during the transpose update");
		for (int c = 0; c < _num_block_cols; ++c) {
			ColumnsMap& transposed_block_row = result_->_block_rows[c];

			ColumnsMap::const_iterator it = block_row.find(c);
			if (it == block_row.end())
				continue;
			else{
				(*result_->_storage[Association(c,r)]) = it->second->transpose();
				transposed_block_row[r] = result_->_storage[Association(c,r)];
			}
		}
	}
}


void SparseBlockMatrix::solveLinearSystem(DenseBlockVector& rhs_vector_,
		DenseBlockVector& result_) const {
	if (result_.size != rhs_vector_.size) {
		result_.reset();
		result_.init(rhs_vector_.size);
	} else {
		result_.clear();
	}

	if(_num_block_rows != _num_block_cols)
		throw std::runtime_error("Error, non squared matrix :(");

	SparseBlockMatrix* L = cholesky();
	SparseBlockMatrix* U = L->transpose();

	DenseBlockVector y;
	L->forwSub(rhs_vector_, y);
	U->backSub(y, result_);
	y.reset();

	delete U;
	delete L;
}


void SparseBlockMatrix::forwSub(DenseBlockVector& rhs_vector_, DenseBlockVector& result_) const {
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


void SparseBlockMatrix::backSub(DenseBlockVector& rhs_vector_, DenseBlockVector& result_) const {
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
			printBlock(i,j);
	}
}

SparseMatrixBlock SparseBlockMatrix::scalarProd(const ColumnsMap& row1_,
		const ColumnsMap& row2_,
		const int max_pos_) const {
	typename ColumnsMap::const_iterator it_1 = row1_.begin();
	typename ColumnsMap::const_iterator it_2 = row2_.begin();

	SparseMatrixBlock result = SparseMatrixBlock::Zero();
	while(it_1 != row1_.end() && it_2 != row2_.end()){
		int col_idx_1 = it_1->first;
		int col_idx_2 = it_2->first;
		if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_){
			return result;
		}
		if(col_idx_1 == col_idx_2){
			SparseMatrixBlock a = (*it_1->second);
			cerr << "a" << endl;
			SparseMatrixBlock b = it_2->second->transpose();
			cerr << "b" << endl;
			result.noalias() += a*b;
//			result += (*it_1->second) * it_2->second->transpose();
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

} /* namespace sparse */

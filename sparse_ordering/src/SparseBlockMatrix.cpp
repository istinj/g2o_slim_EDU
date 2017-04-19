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
}

SparseBlockMatrix::SparseBlockMatrix(const VerticesContainer& vertices_,
		const FactorsMap& factors_){
	_num_block_rows = vertices_.size();
	_num_block_cols = vertices_.size();
	_block_rows.resize(_num_block_rows);

	for (FactorsMap::const_iterator factor = factors_.begin(); factor != factors_.end(); ++factor) {
		setBlock(factor->first.from, factor->first.to, factor->second);
	}
	_has_storage = false;
}

SparseBlockMatrix::SparseBlockMatrix(const std::vector<Vertex>& vertices_,
		const FactorsMap& factors_,
		const std::vector<int> ordering_){
	_num_block_rows = vertices_.size();
	_num_block_cols = vertices_.size();
	_block_rows.resize(_num_block_rows);
	_has_storage = false;
	//! TODO ORDERING
}

//! TODO: how to delete this shit when it owns the memory?
SparseBlockMatrix::~SparseBlockMatrix() {
	//! Add something?
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

	SparseBlockMatrix* chol = new SparseBlockMatrix(_num_block_rows, _num_block_cols, true);
	std::vector<SparseMatrixBlock, Eigen::aligned_allocator<SparseMatrixBlock> > inverse_diag_blocks(_num_block_rows);

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		ColumnsMap& chol_block_row = chol->_block_rows[r];

		if(block_row.empty())
			throw std::runtime_error("Something went wrong during the cholesky");
		int starting_col_idx = block_row.begin()->first;

		for (int c = starting_col_idx; c <= r; ++c) {
			std::shared_ptr<SparseMatrixBlock> accumulator_ptr = std::make_shared<SparseMatrixBlock>();
			std::shared_ptr<SparseMatrixBlock> chol_computed_block = std::make_shared<SparseMatrixBlock>();

			accumulator_ptr->setZero();
			chol_computed_block->setZero();

			ColumnsMap& chol_upper_row = chol->_block_rows[c];
			(*accumulator_ptr) = getBlock(r,c) - (*scalarProd(chol_block_row, chol_upper_row, c-1));

			if (r == c) {
				(*chol_computed_block) = (*accumulator_ptr).llt().matrixL();
				inverse_diag_blocks[r] = (*chol_computed_block).inverse().transpose();
			} else {
				(*chol_computed_block) = (*accumulator_ptr) * inverse_diag_blocks[c];
			}
			chol_block_row[c] = chol_computed_block.get();
		}
	}

	return chol;
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
				std::shared_ptr<SparseMatrixBlock> block = std::make_shared<SparseMatrixBlock>();
				(*block) = it->second->transpose();
				transposed_block_row[r] = block.get();
//				transposed_block_row[r] = it->second;
			}
		}
	}
	return transposed;
}

//! TODO: This shit produces memory leak
SparseBlockMatrix* SparseBlockMatrix::rightMultiplySparseMatrix(const SparseBlockMatrix* other_) const {
	if(other_->_num_block_rows != _num_block_cols)
		throw std::runtime_error("Error, matrices dimensions must agree");
	SparseBlockMatrix* result = new SparseBlockMatrix(_num_block_rows, _num_block_cols, true);
	SparseBlockMatrix* other_transposed = other_->transpose();

	for (int r = 0; r < _num_block_rows; ++r) {
		const ColumnsMap& block_row = _block_rows[r];
		for (int c = 0; c < _num_block_cols; ++c) {
			ColumnsMap& other_transposed_block_row = other_transposed->_block_rows[c];
			shared_ptr<SparseMatrixBlock> result_block = scalarProd(block_row,
					other_transposed_block_row, _num_block_cols);
			result->setBlock(r,c,result_block.get());
		}
	}
	delete other_transposed;
	return result;
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

std::shared_ptr<SparseMatrixBlock> SparseBlockMatrix::scalarProd(const ColumnsMap& row1_,
		const ColumnsMap& row2_,
		const int max_pos_) const{
	ColumnsMap::const_iterator it1 = row1_.begin();
	ColumnsMap::const_iterator it2 = row2_.begin();
	shared_ptr<SparseMatrixBlock> result = make_shared<SparseMatrixBlock>();
	while(it1 != row1_.end() && it2 != row2_.end()){
		int col_idx_1 = it1->first;
		int col_idx_2 = it2->first;
		if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_){
			return result;
		}
		if(col_idx_1 == col_idx_2){
			//! TODO this shit sucks
			(*result) += (*it1->second) * (*it2->second).transpose();
			++it1;
			++it2;
		}
		else if(col_idx_1 > col_idx_2)
			++it2;
		else if(col_idx_1 < col_idx_2)
			++it1;
	}
	return result;
}

SparseMatrixBlock* SparseBlockMatrix::scalarProdPtr(const ColumnsMap& row1_,
		const ColumnsMap& row2_,
		const int max_pos_) const {
	ColumnsMap::const_iterator it1 = row1_.begin();
	ColumnsMap::const_iterator it2 = row2_.begin();
	SparseMatrixBlock* result = new SparseMatrixBlock();
	while(it1 != row1_.end() && it2 != row2_.end()){
		int col_idx_1 = it1->first;
		int col_idx_2 = it2->first;
		if(col_idx_1 > max_pos_ || col_idx_2 > max_pos_){
			return result;
		}
		if(col_idx_1 == col_idx_2){
			//! TODO this shit sucks
			(*result) += (*it1->second) * (*it2->second).transpose();
			++it1;
			++it2;
		}
		else if(col_idx_1 > col_idx_2)
			++it2;
		else if(col_idx_1 < col_idx_2)
			++it1;
	}
	return result;
}

} /* namespace sparse */

using namespace std;

namespace sparse {

template<typename VectorBlockType_>
DenseVector<VectorBlockType_>::DenseVector()
{
   _num_block_rows = 0;
   _total_rows = _num_block_rows * _block_dim;

}

template<typename VectorBlockType_>
DenseVector<VectorBlockType_>::DenseVector(const int num_block_rows_, const int block_dim_){
   _num_block_rows = num_block_rows_;
   _block_dim = block_dim_;
   _total_rows = _num_block_rows * _block_dim;

   _rhs_block_container.resize(num_block_rows_);
}

template<typename VectorBlockType_>
DenseVector<VectorBlockType_>::~DenseVector()
{
   // TODO Auto-generated destructor stub
}

template<typename VectorBlockType_>
void DenseVector<VectorBlockType_>::resize(const int new_rows_){
   _num_block_rows = new_rows_;
   _total_rows = new_rows_ * _block_dim;

   _rhs_block_container.resize(new_rows_);
}

template<typename VectorBlockType_>
void DenseVector<VectorBlockType_>::reset(){
   for(int i = 0; i < _total_rows; ++i){
      _rhs_block_container[i].clear();
   }
   _total_rows = 0;
   _num_block_rows = 0;
}

template<typename VectorBlockType_>
void DenseVector<VectorBlockType_>::setBlock(const int r_, DenseVectorBlock block_){
   if(r_ >= _num_block_rows)
      throw std::runtime_error("set Out of bound");

   if(!block_.isZero()){
      _rhs_block_container[r_] = block_;
   } else {
      _rhs_block_container[r_].setZero();
   }
}

template<typename VectorBlockType_>
VectorBlockType_ DenseVector<VectorBlockType_>::getBlock(const int r_) const {
   if(r_ >= _num_block_rows)
      throw std::runtime_error("get Out of bound");

   return _rhs_block_container[r_];
}



template<typename VectorBlockType_>
void DenseVector<VectorBlockType_>::printBlock(const int r_) const {
   if(r_ >= _num_block_rows)
      throw std::runtime_error("print out of bound");

   cerr << BOLDWHITE << "Block(" << r_ << ")" << ":\n" <<
         CYAN << getBlock(r_) << RESET << endl;
}

template<typename VectorBlockType_>
void DenseVector<VectorBlockType_>::printVector(void) const {
   cerr << BOLDWHITE <<  "Printing vector..." << RESET << endl;
   for (int i = 0; i < numRows(); ++i)
      printBlock(i);
}




} /* namespace sparse */


#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <Eigen/Core>
#include <boost/unordered_map.hpp>

#include "utilities.h"


//! Generic sparse matrix stored in a effient way
namespace sparse {
typedef std::map<int, float> ColumnsMap;
typedef std::vector<ColumnsMap> RowContainer;
class SparseMatrix
{
public:
   SparseMatrix();
   SparseMatrix(int num_rows_, int num_cols_);
   virtual ~SparseMatrix();

   void setElement(const int r_, const int c_, const float value_);
   void resetMatrix(void);
   void resize(const int new_rows_, const int new_cols_);
   void loadFromTXT(const std::string& path_to_file_);

   void evaluateCholeskyStructure(SparseMatrix& cholesky_);
   bool evaluateScalarProdStructure(const ColumnsMap& row_1_,
         const ColumnsMap& row_2_, int max_pos_);

   void cholesky(SparseMatrix& cholesky_);
   float scalarProd(const ColumnsMap& row_1_,
         const ColumnsMap& row_2_, const int max_pos_);

   void printElement(const int r_, const int c_) const;
   bool isNonZeroElement(const int r_, const int c_) const;
   float getElement(const int r_, const int c_) const;

   inline const int numRows(void){return _num_rows;};
   inline const int numCols(void){return _num_cols;};

protected:
   int _num_rows, _num_cols;
   RowContainer _row_container;
};
} // end namespace

#endif

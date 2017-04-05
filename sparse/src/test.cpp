#include <iostream>
#include <Eigen/Core>

#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"
#include "RHSVector.h"

int main(int argc, char const *argv[])
{
   /*
    *
    * A matrix:
   1.3137    1.9933         0         0    1.4421    1.7743
   1.9933    3.1610         0         0    2.3542    2.6394
        0         0    0.2354    0.3563         0         0
        0         0    0.3563    0.5396         0         0
   1.4421    2.3542         0         0    2.1761    2.1368
   1.7743    2.6394         0         0    2.1368    2.6851

    * B RHS vector
    0.5377
    1.8339
   -2.2588
    0.8622
    0.3188
   -1.3077
    *
    *
   /**/


   sparse::SparseBlockMatrix<Eigen::Matrix2d> sp_block_mat(3,3,2);
   sparse::DenseVector<Eigen::Vector2d> rhs_B(3,2);
   Eigen::Matrix2d a;
   Eigen::Vector2d b;

   //! Setup the linear system (retarded fashion)
   a << 1.3137,   1.9933,
         1.9933, 3.1610;
   sp_block_mat.setBlock(0,0,a);

   a << 1.4421, 1.7743,
         2.3542, 2.6394;
   sp_block_mat.setBlock(0,2,a);

   a << 0.2354, 0.3563,
         0.3563, 0.5396;
   sp_block_mat.setBlock(1,1,a);

   a << 1.4421, 2.3542,
         1.7743, 2.6394;
   sp_block_mat.setBlock(2,0,a);

   a << 2.1761, 2.1368,
         2.1368, 2.6851;
   sp_block_mat.setBlock(2,2,a);

   //! RHS
   b << 0.5377,
         1.8339;
   rhs_B.setBlock(0, b);

   b << -2.2588,
         0.8622;
   rhs_B.setBlock(1, b);

   b << 0.3188,
         -1.3077;
   rhs_B.setBlock(2, b);

   sp_block_mat.exportToTxt("../data/cpp_L_export.txt");


   //! Solving linear system Ax = B
//   sparse::DenseVector<Eigen::Vector2d> deltaX(3,2);
//   sp_block_mat.solveLinearSystem(rhs_B, deltaX);

//   cout << endl << BOLDGREEN << "Solution of the system:" << RESET <<  endl;
//   deltaX.printVector();



   return 0;
}


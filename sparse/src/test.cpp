#include <iostream>
#include <Eigen/Core>

#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"
#include "RHSVector.h"

int main(int argc, char const *argv[])
{

	sparse::SparseBlockMatrix<Matrix6f> sp_block_mat(3,3,6);
	sparse::DenseVector<Vector6f> rhs_B(11,6);
	Matrix6f a, a_plus_a;
	Vector6f b;

	//! Setup the linear system (retarded fashion)
	a <<     	163,   118,   135,   107,   190,   144,
			118,   139,   106,   123,   186,   145,
			135,   106,   123,    99,   172,   104,
			107,   123,    99,   155,   202,   117,
			190,   186,   172,   202,   312,   183,
			144,   145,   104,   117,   183,   206;
	a_plus_a <<    326,   236,   270,   214,   380,   288,
			236,   278,   212,   246,   372,   290,
			270,   212,   246,   198,   344,   208,
			214,   246,   198,   310,   404,   234,
			380,   372,   344,   404,   624,   366,
			288,   290,   208,   234,   366,   412;
	for (int i = 0; i < 3; ++i) {
		if(i == 0){
//			a += Matrix6f::Identity() * 1000;
			sp_block_mat.setBlock(i,i,a);
			continue;
		}
		sp_block_mat.setBlock(i,i,a_plus_a);
		if(i > 0)
			sp_block_mat.setBlock(i,i-1,-a);
	}

	sparse::SparseBlockMatrix<Matrix6f> L, U, result;
	sp_block_mat.cholesky(L);
	L.transpose(U);
	L.rightMultiplyMatrix(U, result);
	cout << BOLDGREEN << "A:\n" << RESET << endl;
	sp_block_mat.printMatrix();
	cout << BOLDGREEN << "L*U:\n" << RESET << endl;
	result.printMatrix();



/*
	//! RHS
	b << 3,1,3,1,6,8;
	rhs_B.setBlock(0, b);

	b << 7,1,1,8,10,6;
	rhs_B.setBlock(1, b);

	b << 2,9,4,3,8,1;
	rhs_B.setBlock(2, b);

	b << 1,7,7,6,8,8;
	rhs_B.setBlock(3, b);

	b <<    8,
			3,
			7,
			6,
			4,
			1;
	rhs_B.setBlock(4, b);

	b <<     8,
			4,
			7,
			8,
			2,
			2;
	rhs_B.setBlock(5, b);

	b <<     6,
			5,
			9,
			8,
			8,
			1;
	rhs_B.setBlock(6, b);

	b <<  1,
			1,
			8,
			10,
			7,
			2;
	rhs_B.setBlock(7, b);

	b <<     8,
			2,
			2,
			7,
			4,
			7;
	rhs_B.setBlock(8, b);

	b <<       8,
			6,
			8,
			3,
			8,
			10;
	rhs_B.setBlock(9, b);

	b <<       9,
			1,
			4,
			4,
			7,
			6;
	rhs_B.setBlock(10, b);


	//	//! Test for BW/FW sub
	//	sparse::SparseBlockMatrix<Matrix6f> true_L(11,11,6);
	//	sparse::SparseBlockMatrix<Matrix6f> true_U(11,11,6);
	//	Matrix6f l00, l11;
	//	l00 <<    		   12.7671,         0,         0,        0,         0,         0,
	//			9.2425,    7.3196,         0,        0,         0,         0,
	//			10.5740,    1.1298,    3.1486,        0,         0,         0,
	//			8.3809,    6.2216,    1.0643,   6.7022,         0,         0,
	//			14.8819,    6.6198,    2.2737,   5.0238,    4.0371,         0,
	//			11.2790,    5.5679,   -6.8458,  -0.7286,   -0.6157,    0.0987;
	//	l11 <<     18.0555,         0,         0,         0,         0,         0,
	//			13.0708,   10.3515,         0,         0,         0,         0,
	//			14.9539,    1.5978,    4.4528,         0,         0,         0,
	//			11.8524,    8.7987,    1.5051,    9.4783,         0,         0,
	//			21.0463,    9.3617,    3.2154,    7.1048,    5.7093,         0,
	//			15.9508,    7.8742,   -9.6814,   -1.0304,   -0.8707,    0.1396;
	//
	//	for (int i = 0; i < 11; ++i) {
	//		if(i == 0){
	//			true_L.setBlock(i,i,l00);
	//			continue;
	//		}
	//		true_L.setBlock(i,i,l11);
	//	}

	//	true_L.transpose(true_U);
	//	sparse::DenseVector<Vector6f> Y;
	//	true_L.forwSubstitution(rhs_B, Y);
	//	cout << endl << BOLDGREEN << "Y vector:" << RESET <<  endl;
	//	Y.printVector();
	//	true_U.backSubstitution(Y, deltaX);



	//! Solving linear system Ax = B
//	sparse::DenseVector<Vector6f> deltaX(11,6);
//	sp_block_mat.solveLinearSystem(rhs_B, deltaX);
//
//	cout << endl << BOLDGREEN << "Solution of the system:" << RESET <<  endl;
//	deltaX.printVector();


/**/
	return 0;
}


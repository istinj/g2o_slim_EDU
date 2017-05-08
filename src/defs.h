#pragma once
#include <iostream>
#include <vector>
#include <memory>
#include <ctime>
#include <chrono>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>
#include <suitesparse/cs.h>
#include <suitesparse/amd.h>
#include <suitesparse/cholmod.h>
#include <suitesparse/colamd.h>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#define UNDERLINEDBLACK   "\033[4m\033[30m"      /* Bold Underlined Black */
#define UNDERLINEDRED     "\033[4m\033[31m"      /* Bold Underlined Red */
#define UNDERLINEDGREEN   "\033[4m\033[32m"      /* Bold Underlined Green */
#define UNDERLINEDYELLOW  "\033[4m\033[33m"      /* Bold Underlined Yellow */
#define UNDERLINEDBLUE    "\033[4m\033[34m"      /* Bold Underlined Blue */
#define UNDERLINEDMAGENTA "\033[4m\033[35m"      /* Bold Underlined Magenta */
#define UNDERLINEDCYAN    "\033[4m\033[36m"      /* Bold Underlined Cyan */
#define UNDERLINEDWHITE   "\033[4m\033[37m"      /* Bold Underlined White */

//! For the graph
#define VERTEX_3F "VERTEX_TRACKXYZ"
#define VERTEX_SE3 "VERTEX_SE3:QUAT"
#define EDGE_POSE_POINT "EDGE_SE3_TRACKXYZ"
#define EDGE_SE3 "EDGE_SE3:QUAT"

#define X_DIM 6
#define L_DIM 6

#define PRINT_VAR(x) std::cerr << #x << std::endl << x << std::endl;

//! For the optimizer
#define PLAIN   (uint8_t)0
#define AMD     (uint8_t)1

typedef double Real;
typedef int Counter;
typedef uint8_t SolverType;

typedef std::vector<int> IntVector;
typedef std::pair<int, int> IntPair;

typedef Eigen::Matrix<Real, 2, 2> Matrix2;
typedef Eigen::Matrix<Real, 3, 3> Matrix3;
typedef Eigen::Matrix<Real, 6, 6> Matrix6;
typedef Eigen::Matrix<Real, 3, 6> Matrix3_6;
typedef Eigen::Matrix<Real, 6, 3> Matrix6_3;
typedef Eigen::Matrix<Real, 2, 3> Matrix2_3;
typedef Eigen::Matrix<Real, 2, 6> Matrix2_6;
typedef Eigen::Matrix<Real, 12, 6> Matrix12_6;
typedef Eigen::Matrix<Real, 12, 12> Matrix12;
typedef Eigen::Matrix<Real, 6, 1> Vector6;
typedef Eigen::Matrix<Real, 3, 1> Vector3;
typedef Eigen::Matrix<Real, 2, 1> Vector2;
typedef Eigen::Matrix<Real, 9, 1> Vector9;
typedef Eigen::Matrix<Real, 12, 1> Vector12;
typedef Eigen::Transform<Real,3,Eigen::Isometry> Isometry3;
typedef Eigen::AngleAxis<Real> AngleAxisReal;
typedef Eigen::Quaternion<Real> QuaternionReal;

typedef Vector3 PointMeas;
typedef Isometry3 PoseMeas;
typedef Matrix3 OmegaPoint;
typedef Matrix6 OmegaPose;

typedef Isometry3 Pose;
typedef Vector3 PointXYZ;

typedef Matrix6 SparseMatrixBlock;

typedef Vector6 DenseVectorBlock;
typedef std::map<int, DenseVectorBlock*, std::less<int>,
    Eigen::aligned_allocator<std::pair<int, DenseVectorBlock*> > > IntDenseVectorBlockPtrMap;

typedef std::map<int, SparseMatrixBlock*, std::less<int>,
        Eigen::aligned_allocator<std::pair<const int, SparseMatrixBlock*> > > IntSparseMatrixBlockPtrMap;
typedef std::map<IntPair, SparseMatrixBlock*, std::less<std::pair<int, int> >,
        Eigen::aligned_allocator<std::pair<IntPair, SparseMatrixBlock*> > > IndicesSparseMatrixBlocksPtrMap;






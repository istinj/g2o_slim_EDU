#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Geometry>

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

#define X_DIM 6
#define L_DIM 6

#define PRINT_VAR(x) std::cout << #x << std::endl << x << std::endl;

typedef float real_;

typedef Eigen::Matrix<real_, 2, 2> Matrix2;
typedef Eigen::Matrix<real_, 3, 3> Matrix3;
typedef Eigen::Matrix<real_, 6, 6> Matrix6;
typedef Eigen::Matrix<real_, 3, 6> Matrix3_6;
typedef Eigen::Matrix<real_, 6, 3> Matrix6_3;
typedef Eigen::Matrix<real_, 2, 3> Matrix2_3;
typedef Eigen::Matrix<real_, 2, 6> Matrix2_6;
typedef Eigen::Matrix<real_, 12, 6> Matrix12_6;
typedef Eigen::Matrix<real_, 6, 1> Vector6;
typedef Eigen::Matrix<real_, 3, 1> Vector3;
typedef Eigen::Matrix<real_, 2, 1> Vector2;
typedef Eigen::Matrix<real_, 12, 1> Vector12;
typedef Eigen::Transform<real_,3,Eigen::Isometry> Isometry3;
typedef Eigen::AngleAxis<real_> AngleAxisReal;

typedef Vector3 PointMeas;
typedef Isometry3 PoseMeas;
typedef Matrix3 OmegaPoint;
typedef Matrix6 OmegaPose;

typedef Isometry3 Pose;
typedef Vector3 PointXYZ;

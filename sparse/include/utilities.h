#pragma once
#include <iostream>
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

typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef Eigen::Matrix<float, 6, 6> Matrix6f;

struct Association{
   int x_idx;
   int h_idx;
};

Eigen::Matrix4f v2t(const Vector6f& v);
Eigen::Matrix3f skew(const Eigen::Vector3f& p);

void cholesky(const Matrix6f& input_, Matrix6f& L_);
void loadMatrix(const std::string& name_, Eigen::MatrixXf& data_);

//template<int _Dim>
//void generalCholesky(void){
//	static const int dim = _Dim;
//	Eigen::Matrix<float, dim, dim> A;
//}

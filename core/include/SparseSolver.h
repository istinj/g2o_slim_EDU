/*
 * SparseSolver.h
 *
 *  Created on: 07/mar/2017
 *      Author: istin
 */

#ifndef SPARSESOLVER_H_
#define SPARSESOLVER_H_

#include <iostream>
#include <unordered_map>
#include <vector>
#include <set>

#include <boost/unordered_map.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/Cholesky>

#include "utilities.h"
#include "SparseBlockMatrix.h"
#include "RHSVector.h"
#include "Graph.h"

namespace optimizer {
typedef Eigen::Matrix<float, 6, 6> Matrix6f;
typedef Eigen::Matrix<float, 3, 6> Matrix3_6f;
typedef Eigen::Matrix<float, 6, 3> Matrix6_3f;
typedef Eigen::Matrix<float, 2, 3> Matrix2_3f;
typedef Eigen::Matrix<float, 2, 6> Matrix2_6f;
typedef Eigen::Matrix<float, 12, 6> Matrix12_6f;
typedef Eigen::Matrix<float, 6, 1> Vector6f;
typedef Eigen::Matrix<float, 12, 1> Vector12f;

typedef std::vector<VertexSE3> PosesContainer;
typedef std::vector<VertexXYZ> LandmarkPointsContainer;
typedef std::vector<EdgePosePose> PosePoseEdgeContainer;
typedef std::vector<EdgePosePoint> PosePointEdgeContainer;

typedef std::pair<int, int> HessianIndices;

class SparseSolver {
public:
	SparseSolver();
	SparseSolver(const PosesContainer& robot_poses_,
			const LandmarkPointsContainer& land_points_,
			const PosePoseEdgeContainer& zr_,
			const PosePointEdgeContainer& zl_,
			const float l_, const float epsilon_);
	virtual ~SparseSolver();

	void updateGraph(Graph& graph_);
	void oneStep(void);

private:
	bool linearizePosePoint(float& total_chi_, int& inliers_);
	void linearizePosePose(float& total_chi_, int& inliers_);
	void errorAndJacobianPosePoint(const Pose& xr,
			const PointXYZ& xl,
			const PointMeas& zl,
			Eigen::Vector3f& error,
			Eigen::Matrix3f& Jl,
			Matrix3_6f& Jr);
	void errorAndJacobianPosePose(const Pose& xi,
			const Pose& xj,
			const PoseMeas& zr,
			Vector12f& error,
			Matrix12_6f& Ji,
			Matrix12_6f& Jj);

	inline Eigen::Matrix3f skew(const Eigen::Vector3f& p)
	{
		Eigen::Matrix3f s;
		s <<	0,  -p.z(), p.y(),
				p.z(), 0,  -p.x(),
				-p.y(), p.x(), 0;
		return s;
	}

	inline Pose v2t(const Vector6f& v){
		Pose T = Pose::Identity();
		Eigen::Matrix3f Rx, Ry, Rz;
		Rx = Eigen::AngleAxisf(v(3), Eigen::Vector3f::UnitX());
		Ry = Eigen::AngleAxisf(v(4), Eigen::Vector3f::UnitY());
		Rz = Eigen::AngleAxisf(v(5), Eigen::Vector3f::UnitZ());
		T.linear() = Rx * Ry * Rz;
		T.translation() = v.block<3,1>(0,0);
		return T;
	}

	PosesContainer _robot_poses;
	LandmarkPointsContainer _land_points;

	PosePoseEdgeContainer _Zr;
	PosePointEdgeContainer _Zl;

	//! TODO Remember to clean-up everything in the destructor (or at the end of the iteration)
	//! TODO Same for the RHSVector;
	sparse::SparseBlockMatrix<Matrix6f>* _pose_pose_Hessian;
	sparse::DenseVector<Vector6f>* _pose_pose_B;

	float _lambda = 0.0;
	float _threshold = 0.0;


public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace optimizer */

#endif /* SPARSESOLVER_H_ */

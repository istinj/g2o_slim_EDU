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

#include "defs.h"
#include "SparseBlockMatrix.h"
#include "RHSVector.h"
#include "Graph.h"

namespace optimizer {
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
			const real_ l_, const real_ epsilon_);
	virtual ~SparseSolver();

	void updateGraph(Graph& graph_);
	void oneStep(void);

private:
	bool linearizePosePoint(real_& total_chi_, int& inliers_);
	void linearizePosePose(real_& total_chi_, int& inliers_);
	void errorAndJacobianPosePoint(const Pose& xr,
			const PointXYZ& xl,
			const PointMeas& zl,
			Vector3& error,
			Matrix3& Jl,
			Matrix3_6& Jr);
	void errorAndJacobianPosePose(const Pose& xi,
			const Pose& xj,
			const PoseMeas& zr,
			Vector12& error,
			Matrix12_6& Ji,
			Matrix12_6& Jj);

	inline Matrix3 skew(const Vector3& p)
	{
		Matrix3 s;
		s <<	0,  -p.z(), p.y(),
				p.z(), 0,  -p.x(),
				-p.y(), p.x(), 0;
		return s;
	}

	inline Pose v2t(const Vector6& v){
		Pose T = Pose::Identity();
		Matrix3 Rx, Ry, Rz;
		Rx = AngleAxisReal(v(3), Vector3::UnitX());
		Ry = AngleAxisReal(v(4), Vector3::UnitY());
		Rz = AngleAxisReal(v(5), Vector3::UnitZ());
		T.linear() = Rx * Ry * Rz;
		T.translation() = v.block<3,1>(0,0);
		return T;
	}

	PosesContainer _robot_poses;
	LandmarkPointsContainer _land_points;

	PosePoseEdgeContainer _Zr;
	PosePointEdgeContainer _Zl;

	sparse::SparseBlockMatrix<Matrix6>* _pose_pose_Hessian;
	sparse::DenseVector<Vector6>* _pose_pose_B;

	real_ _lambda = 0.0;
	real_ _threshold = 0.0;

	Matrix3 _Rx0, _Ry0, _Rz0;


public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace optimizer */

#endif /* SPARSESOLVER_H_ */

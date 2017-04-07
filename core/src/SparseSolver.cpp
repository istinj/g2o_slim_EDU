/*
 * SparseSolver.cpp
 *
 *  Created on: 07/mar/2017
 *      Author: istin
 */

#include "SparseSolver.h"

using namespace std;
using namespace Eigen;

namespace optimizer {

SparseSolver::SparseSolver() {
	_pose_pose_Hessian = new sparse::SparseBlockMatrix<Matrix6f>(1,1,6);
	_pose_pose_B = new sparse::DenseVector<Vector6f>(1,6);
}

SparseSolver::~SparseSolver() {
	// TODO Auto-generated destructor stub
	delete _pose_pose_B;
	delete _pose_pose_Hessian;
}

SparseSolver::SparseSolver(const PosesContainer& robot_poses_,
		const LandmarkPointsContainer& land_points_,
		const PosePoseEdgeContainer& zr_,
		const PosePointEdgeContainer& zl_,
		const float l_, const float epsilon_){
	_robot_poses = robot_poses_;
	_land_points = land_points_;
	_Zr = zr_;
	_Zl = zl_;
	_lambda = l_;
	_threshold = epsilon_;

	_pose_pose_Hessian = new sparse::SparseBlockMatrix<Matrix6f>(robot_poses_.size(),
			robot_poses_.size(),6);
	_pose_pose_B = new sparse::DenseVector<Vector6f>(robot_poses_.size(),6);
}

bool SparseSolver::linearizePosePoint(float& total_chi_, int& inliers_){
	Matrix3f Jl = Matrix3f::Zero();
	Matrix3_6f Jr = Matrix3_6f::Zero();
	Vector3f e = Vector3f::Zero();

	total_chi_= 0.0;
	inliers_ = 0;

	for(PosePointEdgeContainer::iterator it = _Zl.begin(); it != _Zl.end(); ++it){
		std::pair<int, int> curr_association = it->association();

		PosesContainer::iterator pose_iter = std::find(_robot_poses.begin(),
				_robot_poses.end(), curr_association.first);
		if(pose_iter == _robot_poses.end()){
			cerr << "Error, bad data association" << endl;
			return false;
		}
		LandmarkPointsContainer::iterator land_iter = std::find(_land_points.begin(),
				_land_points.end(), curr_association.second);
		if(land_iter == _land_points.end()){
			cerr << "Error, bad data association" << endl;
			return false;
		}
		errorAndJacobianPosePoint(pose_iter->data(),
				land_iter->data(),
				it->data(),
				e, Jl, Jr);

		float chi = e.transpose() * e;
		if(chi > _threshold){
			e *= sqrt(_threshold/chi);
			chi = _threshold;
		} else {
			inliers_++;
		}
		total_chi_ += chi;


		//! TODO H?
		int pose_hessian_idx = pose_iter->index();
		int land_hessian_idx = land_iter->index();

		Matrix6f h_pp_data = Jr.transpose() * Jr;
		pair<int, int> hessian_indices = make_pair(pose_hessian_idx,
				pose_hessian_idx);

		//! This block is useless
		Matrix6_3f h_pl_data = Jr.transpose() * Jl;
		hessian_indices = make_pair(pose_hessian_idx,
				land_hessian_idx);

		Matrix3_6f h_lp_data = Jl.transpose() * Jr;
		hessian_indices = make_pair(land_hessian_idx,
				pose_hessian_idx);

		Matrix3f h_ll_data = Jl.transpose() * Jl;
		hessian_indices = make_pair(land_hessian_idx,
				land_hessian_idx);

		//! TODO B?
/*
		b_block_pose.data = Jr.transpose() * e;
		b_block_pose.blockIndex = pose_hessian_idx;

		b_block_land.data = Jl.transpose() * e;
		b_block_land.blockIndex = pose_hessian_idx;
/**/
	}
	return true;
}

bool SparseSolver::linearizePosePose(float& total_chi_, int& inliers_) {
	total_chi_ = 0.0;
	inliers_ = 0;

	Matrix12_6f Ji = Matrix12_6f::Zero();
	Matrix12_6f Jj = Matrix12_6f::Zero();
	Vector12f e = Vector12f::Zero();

	Matrix<float, 12, 12> Omega = Matrix<float, 12, 12>::Identity();
	Omega.block<9,9>(0,0) *= 1000.0;

	//! For each Pose-Pose edge
	for (PosePoseEdgeContainer::iterator it = _Zr.begin(); it != _Zr.end();	++it) {
		std::pair<int, int> curr_association = it->association();

		//! Retrieve the Pose Vertices involved in the current edge
		PosesContainer::iterator pose_i_iter = std::find(_robot_poses.begin(),
				_robot_poses.end(), curr_association.first);
		if(pose_i_iter == _robot_poses.end()){
			cerr << "Error, bad data association" << endl;
			return false;
		}

		PosesContainer::iterator pose_j_iter = std::find(_robot_poses.begin(),
				_robot_poses.end(), curr_association.second);
		if(pose_j_iter == _robot_poses.end()){
			cerr << "Error, bad data association" << endl;
			return false;
		}

		//! Evaluate error and jacobian
		errorAndJacobianPosePose(pose_i_iter->data(),
				pose_j_iter->data(),
				it->data(),
				e, Ji, Jj);

		//! Evaluate statistics
		float chi = e.transpose() * Omega * e;
		if (chi > _threshold) {
			Omega *= sqrt(_threshold / chi);
			chi = _threshold;
		} else {
			inliers_++;
		}
		total_chi_ += chi;


		//! Start building the Hessian
		//! Each meas introduces 3 blocks in the Hessian Matrix, H_ii, H_ji, H_jj,
		//! Building only the lower triangular part of H;
		//! Block indices
		int i_hessian_idx = pose_i_iter->index();
		int j_hessian_idx = pose_j_iter->index();

		//! First block
		Matrix6f h_ii_block = Ji.transpose() * Omega * Ji;
		HessianIndices h_ii_indices = make_pair(i_hessian_idx,
				i_hessian_idx);
		Matrix6f prev_hii_block = _pose_pose_Hessian->getBlock(h_ii_indices.first,
				h_ii_indices.second);
		h_ii_block.noalias() += prev_hii_block;
		_pose_pose_Hessian->setBlock(h_ii_indices.first,
				h_ii_indices.second, h_ii_block);

		//! Second block
		Matrix6f h_ji_block = Jj.transpose() * Omega * Ji;
		HessianIndices h_ji_indices = make_pair(j_hessian_idx,
				i_hessian_idx);
		Matrix6f prev_h_ji_block = _pose_pose_Hessian->getBlock(h_ji_indices.first,
				h_ji_indices.second);
		h_ji_block.noalias() += prev_h_ji_block;
		_pose_pose_Hessian->setBlock(h_ji_indices.first,
				h_ji_indices.second, h_ji_block);

		//! Second block transposed; this block is useless
		Matrix6f h_ij_block = Ji.transpose() * Omega * Jj;
		HessianIndices h_ij_indices = make_pair(i_hessian_idx,
				j_hessian_idx);
		Matrix6f prev_h_ij_block = _pose_pose_Hessian->getBlock(h_ij_indices.first,
				h_ij_indices.second);
		h_ij_block.noalias() += prev_h_ij_block;
		_pose_pose_Hessian->setBlock(h_ij_indices.first,
				h_ij_indices.second, h_ij_block);

		//! Third block
		Matrix6f h_jj_block = Jj.transpose() * Omega * Jj;
		HessianIndices h_jj_indices = make_pair(j_hessian_idx,
				j_hessian_idx);
		Matrix6f prev_hjj_block = _pose_pose_Hessian->getBlock(h_jj_indices.first,
				h_jj_indices.second);
		h_jj_block.noalias() += prev_hjj_block;
		_pose_pose_Hessian->setBlock(h_jj_indices.first,
				h_jj_indices.second, h_jj_block);

		//! Building the RHS Vector
		Vector6f b_i = Ji.transpose() * Omega * e;
		Vector6f prev_b_i = _pose_pose_B->getBlock(i_hessian_idx);
		b_i.noalias() += prev_b_i;
		_pose_pose_B->setBlock(i_hessian_idx, b_i);

		Vector6f b_j = Jj.transpose() * Omega * e;
		Vector6f prev_b_j = _pose_pose_B->getBlock(j_hessian_idx);
		b_j.noalias() += prev_b_j;
		_pose_pose_B->setBlock(j_hessian_idx, b_j);
		//! Seems to work fine
	}
	return true;
}
void SparseSolver::errorAndJacobianPosePoint(const Pose& xr,
		const PointXYZ& xl,
		const PointMeas& zl,
		Eigen::Vector3f& error,
		Eigen::Matrix3f& Jl,
		Matrix3_6f& Jr){
	Vector3f h_x = xr.linear() * xl + xr.translation();

	error = h_x - zl;

	Jl = xr.linear();
	Jr.block<3,3>(0,0).setIdentity();
	Jr.block<3,3>(0,3) = -skew(h_x);
}

void SparseSolver::errorAndJacobianPosePose(const Pose& xi,
		const Pose& xj,
		const PoseMeas& zr,
		Vector12f& error,
		Matrix12_6f& Ji,
		Matrix12_6f& Jj){

	//! TODO move to initialization of the solver
	Matrix3f Rx0, Ry0, Rz0;
	Rx0 << 0,0,0,  0,0,-1,  0,1,0;
	Ry0 << 0,0,1,  0,0,0,   -1,0,0;
	Rz0 << 0,-1,0, 1,0,0,   0,0,0;

	Matrix3f Ri = xi.linear();
	Matrix3f Rj = xj.linear();
	Vector3f ti = xi.translation();
	Vector3f tj = xj.translation();
	Vector3f t_ij = tj-ti;

	Matrix3f dR_x = Ri.transpose() * Rx0 * Rj;
	Matrix3f dR_y = Ri.transpose() * Ry0 * Rj;
	Matrix3f dR_z = Ri.transpose() * Rz0 * Rj;

	Matrix<float, 9, 1> dr_x_flattened, dr_y_flattened, dr_z_flattened;
	dr_x_flattened << dR_x.col(0), dR_x.col(1), dR_x.col(2);
	dr_y_flattened << dR_y.col(0), dR_y.col(1), dR_y.col(2);
	dr_z_flattened << dR_z.col(0), dR_z.col(1), dR_z.col(2);

	//! Fill Jj
	Jj.block<9,1>(0,3) = dr_x_flattened;
	Jj.block<9,1>(0,4) = dr_y_flattened;
	Jj.block<9,1>(0,5) = dr_z_flattened;
	Jj.block<3,3>(9,0) = Ri.transpose();
	Jj.block<3,3>(9,3) = -Ri.transpose() * skew(tj);

	Ji = -Jj;

	Pose h_x = Pose::Identity();
	h_x.linear() = Ri.transpose() * Rj;
	h_x.translation() = Ri.transpose() * t_ij;

	//! Compose e
	Eigen::Isometry3f temp_e;
	temp_e.matrix() = h_x.matrix() - zr.matrix();

	error.block<3,1>(0,0) = temp_e.matrix().block<3,1>(0,0);
	error.block<3,1>(3,0) = temp_e.matrix().block<3,1>(0,1);
	error.block<3,1>(6,0) = temp_e.matrix().block<3,1>(0,2);
	error.block<3,1>(9,0) = temp_e.matrix().block<3,1>(0,3);
}

void SparseSolver::updateGraph(Graph& graph_){
	cerr << BOLDYELLOW << "Updating the graph..." << RESET << endl;
	graph_.updateVerticesSE3(_robot_poses);
//	graph_.updateVerticesXYZ(_land_points);
	cerr << BOLDGREEN << "Done!" << RESET << endl;
}


void SparseSolver::oneStep(void){
	float step_chi;
	int step_inliers;

	if(linearizePosePose(step_chi,step_inliers)){
		cerr << GREEN << "inliers odom = " << step_inliers << "\t" << "chi odom = " << step_chi << RESET << endl;

		//! Solving the Linear System Hx = B. Since it is under-determined (the solution is up-to a
		//! rigid transformation), it is necessary to introduce a bias on the first element, to fix the
		//! first point. Moreover, since we want that the starting pose remains the same, it is necessary
		//! to set to 0 the first dX block.
		sparse::DenseVector<Vector6f> dX_pose_pose;

		Matrix6f temp = _pose_pose_Hessian->getBlock(0,0);
		temp += Matrix6f::Identity()*100000.0; //! Bias
		_pose_pose_Hessian->setBlock(0,0,temp);

		_pose_pose_Hessian->solveLinearSystem((*_pose_pose_B), dX_pose_pose);

		dX_pose_pose.setBlock(0, Vector6f::Zero()); //! Fix the starting pose;

		//! Apply the dX to the state.
		for (int i = 0; i < dX_pose_pose.numRows(); ++i) {
			Pose new_pose = v2t(dX_pose_pose.getBlock(i)) * _robot_poses[i].data();
			_robot_poses[i].setData(new_pose);
		}
	} else {
		throw std::runtime_error("Linearize Pose-Pose Failure");
	}

	//! TODO CLEAN-UP EVERYTHING
	return; //placeholder
}

int SparseSolver::getPoseMatrixIndex(int curr_pose_idx){
	if(curr_pose_idx > _robot_poses.size() - 1){
		cerr << "Exceeding index\nExit" << endl;
		return -1;
	}
	return curr_pose_idx * X_DIM;
}

int SparseSolver::getLandMatrixIndex(int curr_land_idx){
	if(curr_land_idx > _land_points.size() - 1){
		cerr << "Exceeding index\nExit" << endl;
		return -1;
	}
	return _robot_poses.size() * X_DIM + curr_land_idx * L_DIM;
}

} /* namespace optimizer */

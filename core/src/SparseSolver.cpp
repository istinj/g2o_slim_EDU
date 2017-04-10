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
	_pose_pose_Hessian = new sparse::SparseBlockMatrix<Matrix6>(1,1,6);
	_pose_pose_B = new sparse::DenseVector<Vector6>(1,6);
}

SparseSolver::~SparseSolver() {
	// TODO Auto-generated destructor stub
}

SparseSolver::SparseSolver(const PosesContainer& robot_poses_,
		const LandmarkPointsContainer& land_points_,
		const PosePoseEdgeContainer& zr_,
		const PosePointEdgeContainer& zl_,
		const real_ l_, const real_ epsilon_){
	_robot_poses = robot_poses_;
	_land_points = land_points_;
	_Zr = zr_;
	_Zl = zl_;
	_lambda = l_;
	_threshold = epsilon_;
}

bool SparseSolver::linearizePosePoint(real_& total_chi_, int& inliers_){
	Matrix3 Jl = Matrix3::Zero();
	Matrix3_6 Jr = Matrix3_6::Zero();
	Vector3 e = Vector3::Zero();

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

		real_ chi = e.transpose() * e;
		if(chi > _threshold){
			e *= sqrt(_threshold/chi);
			chi = _threshold;
		} else {
			inliers_++;
		}
		total_chi_ += chi;


		//! Hessain
		int pose_hessian_idx = pose_iter->index();
		int land_hessian_idx = land_iter->index();

		Matrix6 h_pp_data = Jr.transpose() * Jr;
		pair<int, int> hessian_indices = make_pair(pose_hessian_idx,
				pose_hessian_idx);

		//! This block is useless
		Matrix6_3 h_pl_data = Jr.transpose() * Jl;
		hessian_indices = make_pair(pose_hessian_idx,
				land_hessian_idx);

		Matrix3_6 h_lp_data = Jl.transpose() * Jr;
		hessian_indices = make_pair(land_hessian_idx,
				pose_hessian_idx);

		Matrix3 h_ll_data = Jl.transpose() * Jl;
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

void SparseSolver::linearizePosePose(real_& total_chi_, int& inliers_) {
	total_chi_ = 0.0;
	inliers_ = 0;

	Matrix12_6 Ji = Matrix12_6::Zero();
	Matrix12_6 Jj = Matrix12_6::Zero();
	Vector12 e = Vector12::Zero();

	Matrix<real_, 12, 12> Omega = Matrix<real_, 12, 12>::Identity();
	Omega.block<9,9>(0,0) *= 1000.0;

	//! For each Pose-Pose edge
	for (PosePoseEdgeContainer::iterator it = _Zr.begin(); it != _Zr.end();	++it) {
		std::pair<int, int> curr_association = it->association();

		//! Retrieve the Pose Vertices involved in the current edge
		PosesContainer::iterator pose_i_iter = std::find(_robot_poses.begin(),
				_robot_poses.end(), curr_association.first);
		if(pose_i_iter == _robot_poses.end())
			throw std::runtime_error("Error, bad data association");

		PosesContainer::iterator pose_j_iter = std::find(_robot_poses.begin(),
				_robot_poses.end(), curr_association.second);
		if(pose_j_iter == _robot_poses.end())
			throw std::runtime_error("Error, bad data association");

		//! Evaluate error and jacobian
		errorAndJacobianPosePose(pose_i_iter->data(),
				pose_j_iter->data(),
				it->data(),
				e, Ji, Jj);

		//! Evaluate statistics
		real_ chi = e.transpose() * Omega * e;
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
		Matrix6 h_ii_block = Ji.transpose() * Omega * Ji;
		HessianIndices h_ii_indices = make_pair(i_hessian_idx,
				i_hessian_idx);
		Matrix6 prev_hii_block = _pose_pose_Hessian->getBlock(h_ii_indices.first,
				h_ii_indices.second);
		h_ii_block.noalias() += prev_hii_block;
		_pose_pose_Hessian->setBlock(h_ii_indices.first,
				h_ii_indices.second, h_ii_block);

		//! Second block
		Matrix6 h_ji_block = Jj.transpose() * Omega * Ji;
		HessianIndices h_ji_indices = make_pair(j_hessian_idx,
				i_hessian_idx);
		Matrix6 prev_h_ji_block = _pose_pose_Hessian->getBlock(h_ji_indices.first,
				h_ji_indices.second);
		h_ji_block.noalias() += prev_h_ji_block;
		_pose_pose_Hessian->setBlock(h_ji_indices.first,
				h_ji_indices.second, h_ji_block);

		//! Second block transposed; this block is useless
		Matrix6 h_ij_block = Ji.transpose() * Omega * Jj;
		HessianIndices h_ij_indices = make_pair(i_hessian_idx,
				j_hessian_idx);
		Matrix6 prev_h_ij_block = _pose_pose_Hessian->getBlock(h_ij_indices.first,
				h_ij_indices.second);
		h_ij_block.noalias() += prev_h_ij_block;
		_pose_pose_Hessian->setBlock(h_ij_indices.first,
				h_ij_indices.second, h_ij_block);

		//! Third block
		Matrix6 h_jj_block = Jj.transpose() * Omega * Jj;
		HessianIndices h_jj_indices = make_pair(j_hessian_idx,
				j_hessian_idx);
		Matrix6 prev_hjj_block = _pose_pose_Hessian->getBlock(h_jj_indices.first,
				h_jj_indices.second);
		h_jj_block.noalias() += prev_hjj_block;
		_pose_pose_Hessian->setBlock(h_jj_indices.first,
				h_jj_indices.second, h_jj_block);

		//! Building the RHS Vector
		Vector6 b_i = Ji.transpose() * Omega * e;
		Vector6 prev_b_i = _pose_pose_B->getBlock(i_hessian_idx);
		b_i.noalias() += prev_b_i;
		_pose_pose_B->setBlock(i_hessian_idx, b_i);

		Vector6 b_j = Jj.transpose() * Omega * e;
		Vector6 prev_b_j = _pose_pose_B->getBlock(j_hessian_idx);
		b_j.noalias() += prev_b_j;
		_pose_pose_B->setBlock(j_hessian_idx, b_j);
		//! Seems to work fine
	}
}

void SparseSolver::errorAndJacobianPosePoint(const Pose& xr,
		const PointXYZ& xl,
		const PointMeas& zl,
		Vector3& error,
		Matrix3& Jl,
		Matrix3_6& Jr){
//	Vector3f h_x = xr.linear() * xl + xr.translation();
	Pose inv_xr = xr.inverse();
	Vector3 h_x = inv_xr.linear() * xl + inv_xr.translation();

	error = h_x - zl;

//	Jl = xr.linear();
//	Jr.block<3,3>(0,0).setIdentity();
//	Jr.block<3,3>(0,3) = -skew(h_x);

	Jl = inv_xr.linear();
	Jr.block<3,3>(0,0) = -inv_xr.linear();
	Jr.block<3,3>(0,3) = inv_xr.linear() * skew(h_x);
}

void SparseSolver::errorAndJacobianPosePose(const Pose& xi,
		const Pose& xj,
		const PoseMeas& zr,
		Vector12& error,
		Matrix12_6& Ji,
		Matrix12_6& Jj){

	_Rx0 << 0,0,0,  0,0,-1,  0,1,0;
	_Ry0 << 0,0,1,  0,0,0,   -1,0,0;
	_Rz0 << 0,-1,0, 1,0,0,   0,0,0;

	Matrix3 Ri = xi.linear();
	Matrix3 Rj = xj.linear();
	Vector3 ti = xi.translation();
	Vector3 tj = xj.translation();
	Vector3 t_ij = tj-ti;

	Matrix3 dR_x = Ri.transpose() * _Rx0 * Rj;
	Matrix3 dR_y = Ri.transpose() * _Ry0 * Rj;
	Matrix3 dR_z = Ri.transpose() * _Rz0 * Rj;

	Matrix<real_, 9, 1> dr_x_flattened, dr_y_flattened, dr_z_flattened;
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
	Isometry3 temp_e = Isometry3::Identity();
	temp_e.matrix() = h_x.matrix() - zr.matrix();

	error.setZero();
	error.block<3,1>(0,0) = temp_e.matrix().block<3,1>(0,0);
	error.block<3,1>(3,0) = temp_e.matrix().block<3,1>(0,1);
	error.block<3,1>(6,0) = temp_e.matrix().block<3,1>(0,2);
	error.block<3,1>(9,0) = temp_e.matrix().block<3,1>(0,3);
}

//! This fuction updates the graph after the optimization of
//! poses and landmarks. For now it only updates the poses, but
//! there will be a refactoring.
void SparseSolver::updateGraph(Graph& graph_){
	cerr << endl << BOLDYELLOW << "Updating the graph..." << RESET << endl;
	graph_.updateVerticesSE3(_robot_poses);
	cerr << BOLDGREEN << "Done!" << RESET << endl;
}


void SparseSolver::oneStep(void){
	_pose_pose_Hessian = new sparse::SparseBlockMatrix<Matrix6>(_robot_poses.size(),
			_robot_poses.size(),6);
	_pose_pose_B = new sparse::DenseVector<Vector6>(_robot_poses.size(),6);

	real_ step_chi;
	int step_inliers;
	cerr << BOLDYELLOW <<  "One step with damping = " << BOLDRED << _lambda
			<< BOLDYELLOW << " and kernel threshold = " << BOLDRED << _threshold << RESET << endl;
	linearizePosePose(step_chi,step_inliers);
	cerr << BOLDWHITE << "inliers pose-pose = " << BOLDBLUE << step_inliers << "\t" << BOLDWHITE
			<< "chi pose-pose = " << BOLDBLUE<< step_chi << RESET << endl;

	//! Solving the Linear System Hx = B. Since it is under-determined (the solution is up-to a
	//! rigid transformation), it is necessary to introduce a bias on the first element, to fix the
	//! first point. Moreover, since we want that the starting pose remains the same, it is necessary
	//! to set to 0 the first dX block.
	sparse::DenseVector<Vector6> dX_pose_pose;

	Matrix6 temp = _pose_pose_Hessian->getBlock(0,0);
	temp += Matrix6::Identity()*1000000.0; //! Bias
	_pose_pose_Hessian->setBlock(0,0,temp);


	//! Add damping
	if(_lambda != 0)
		_pose_pose_Hessian->addDamping(_lambda);

	//! Actaully solve the linear system
	_pose_pose_Hessian->solveLinearSystem((*_pose_pose_B), dX_pose_pose);

	dX_pose_pose.setBlock(0, Vector6::Zero()); //! Fix the starting pose;

	//! Apply the dX to the state.
	for (int i = 0; i < dX_pose_pose.numRows(); ++i) {
		Pose new_pose = Pose::Identity();
		new_pose = v2t(-dX_pose_pose.getBlock(i)) * _robot_poses[i].data();
		_robot_poses[i].setData(new_pose);
	}
	//! TODO CLEAN-UP EVERYTHING?
	delete _pose_pose_B;
	delete _pose_pose_Hessian;
}

} /* namespace optimizer */

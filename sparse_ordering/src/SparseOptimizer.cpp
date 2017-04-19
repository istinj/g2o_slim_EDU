/*
 * SparseOptimizer.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <SparseOptimizer.h>

using namespace std;

namespace sparse {

SparseOptimizer::SparseOptimizer() {
	// TODO Auto-generated constructor stub

}

SparseOptimizer::~SparseOptimizer() {
	for (FactorsMap::iterator it = _factors.begin(); it != _factors.end(); ++it){
		delete it->second;
	}

	_B.reset();
	_dX.reset();
}

void SparseOptimizer::init(const VerticesContainer& vertices_,
		const EdgesContainer& edges_){
	_vertices = vertices_;
	_edges = edges_;

	//! Initialize the Factors' Storage
	FactorsMap::iterator f_it;
	for (int i = 0; i < _edges.size(); ++i) {
		VerticesContainer::const_iterator pose_i = std::find(_vertices.begin(),
				_vertices.end(), _edges[i].association().first);
		VerticesContainer::const_iterator pose_j = std::find(_vertices.begin(),
				_vertices.end(), _edges[i].association().second);
		int from = pose_i->index();
		int to = pose_j->index();

		//! Storage for the Hessian
		f_it = _factors.find(Factor(from, from));
		if(f_it == _factors.end()){
			_factors[Factor(from, from)] = new SparseMatrixBlock();
			_factors[Factor(from, from)]->setZero();
		}

		f_it = _factors.find(Factor(to, from));
		if(f_it == _factors.end()){
			_factors[Factor(to, from)] = new SparseMatrixBlock();
			_factors[Factor(to, from)]->setZero();
		}

		f_it = _factors.find(Factor(to, to));
		if(f_it == _factors.end()){
			_factors[Factor(to, to)] = new SparseMatrixBlock();
			_factors[Factor(to, to)]->setZero();
		}

	}

	//! Storage for the rhs vector
	_B.init(_vertices.size());

	//! Storage for the increment vector
	_dX.init(_vertices.size());

	//! TODO	Storage for U and L
	return;
}

void SparseOptimizer::oneStep(void){
	real_ step_chi = 0.0;
	int step_inliers = 0;

	linearizeFactor(step_chi, step_inliers);
	cerr << BOLDWHITE << "inliers pose-pose = " << BOLDBLUE << step_inliers << "\t" << BOLDWHITE
			<< "chi pose-pose = " << BOLDBLUE<< step_chi << RESET << endl;

	//! TODO Build and solve the linear system HdX = B;
/*
	hessian = SparseBlockMatrix(_factors, _ordering);
	SparseMatrixBlock& h_00 = _hessian(0,0) + (SparseMatrixBlock::Identity() * 1000.0);

	hessian.solveLinearSystem(_B, dX);
	updateVertices(_dX);

	clearFactors();
	_B.clear();
/**/
	SparseMatrixBlock& h_00 = (*_factors[Factor(0,0)]);
	h_00 += SparseMatrixBlock::Identity() * 1000000.0; //! Not good but ok for now

	SparseBlockMatrix* hessian = new SparseBlockMatrix(_vertices, _factors); //! OK
	hessian->solveLinearSystem(_B, _dX);

	for (int i = 0; i < _dX.size; ++i) {
		Pose new_pose = Pose::Identity();
		new_pose = v2t(-(*_dX.blocks[i])) * _vertices[i].data();
		_vertices[i].setData(new_pose);
	}

	_dX.clear();
	_B.clear();
	delete hessian; //! Does not delete the matrices pull (just 72kB);
	for (FactorsMap::iterator it = _factors.begin(); it != _factors.end(); ++it){
		it->second->setZero();
	}
}
void SparseOptimizer::linearizeFactor(real_& total_chi_, int& inliers_){
	total_chi_ = 0.0;
	inliers_ = 0;

	Matrix12 omega = Matrix12::Identity();
	omega.block<9,9>(0,0) *= 1000.0;

	Matrix12_6 Ji = Matrix12_6::Zero();
	Matrix12_6 Jj = Matrix12_6::Zero();
	Vector12 e = Vector12::Zero();

	//! For each edge, compute - numerically - the relative factor
	for (int i = 0; i < _edges.size(); ++i) {
		const Edge& edge = _edges[i];
		VerticesContainer::const_iterator pose_i = std::find(_vertices.begin(),
				_vertices.end(), edge.association().first);
		VerticesContainer::const_iterator pose_j = std::find(_vertices.begin(),
				_vertices.end(), edge.association().second);

		int from = pose_i->index();
		int to = pose_j->index();

		//! Compute error and jacobian
		errorAndJacobian(pose_i->data(), pose_j->data(), _edges[i].data(),
				e, Ji, Jj);

		//! Robust kernel
		real_ chi = e.transpose() * omega * e;
		if(chi > _kernel_threshold) {
			omega *= sqrtf(_kernel_threshold / chi);
			chi = _kernel_threshold;
		} else {
			inliers_++;
		}
		total_chi_ += chi;

		//! TODO THIS SHIT SUCKS
		//! Compute the factor contribution to the Hessian
		SparseMatrixBlock& H_ii = (*_factors[Factor(from, from)]);
		SparseMatrixBlock& H_ji = (*_factors[Factor(to, from)]);
		SparseMatrixBlock& H_jj = (*_factors[Factor(to, to)]);
		H_ii.noalias() += Ji.transpose() * omega * Ji;
		H_ji.noalias() += Jj.transpose() * omega * Ji;
		H_jj.noalias() += Jj.transpose() * omega * Jj;

		//! Compute the factor contribution to the rhs vector
		DenseVectorBlock& b_i = (*_B.blocks[from]);
		DenseVectorBlock& b_j = (*_B.blocks[to]);
		b_i.noalias() +=  Ji.transpose() * omega * e;
		b_j.noalias() +=  Jj.transpose() * omega * e;
	}
}

void SparseOptimizer::errorAndJacobian(const Pose& xi, const Pose& xj,const PoseMeas& zr,
		Vector12& error, Matrix12_6& Ji, Matrix12_6& Jj){
	Matrix3 Rx0, Ry0, Rz0;
	Rx0 << 0,0,0,  0,0,-1,  0,1,0;
	Ry0 << 0,0,1,  0,0,0,   -1,0,0;
	Rz0 << 0,-1,0, 1,0,0,   0,0,0;

	Matrix3 Ri = xi.linear();
	Matrix3 Rj = xj.linear();
	Vector3 ti = xi.translation();
	Vector3 tj = xj.translation();
	Vector3 t_ij = tj-ti;

	Matrix3 dR_x = Ri.transpose() * Rx0 * Rj;
	Matrix3 dR_y = Ri.transpose() * Ry0 * Rj;
	Matrix3 dR_z = Ri.transpose() * Rz0 * Rj;

	Vector9 dr_x_flattened, dr_y_flattened, dr_z_flattened;
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

} /* namespace sparse */

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
	for (BlocksMap::iterator it = _blocks_pull.begin(); it != _blocks_pull.end(); ++it){
		delete it->second;
	}

	delete _H;
	delete _L;
	delete _U;

	_B.reset();
	_Y.reset();
	_dX.reset();
}

void SparseOptimizer::init(const VerticesContainer& vertices_,
		const EdgesContainer& edges_){
	_vertices = vertices_;
	_edges = edges_;

	//! Initialize the Factors' Storage
	HessianBlocksMap::iterator f_it;
	for (int i = 0; i < _edges.size(); ++i) {
		VerticesContainer::const_iterator pose_i = std::find(_vertices.begin(),
				_vertices.end(), _edges[i].association().first);
		VerticesContainer::const_iterator pose_j = std::find(_vertices.begin(),
				_vertices.end(), _edges[i].association().second);
		int from = pose_i->index();
		int to = pose_j->index();

		cerr << "from: " << from << "\tto: " << to << endl;

		//! Storage for the Hessian
		f_it = _blocks_pull.find(Association(from, from));
		if(f_it == _blocks_pull.end()){
			_blocks_pull[Association(from, from)] = new SparseMatrixBlock();
			_blocks_pull[Association(from, from)]->setIdentity();
			cerr << BOLDGREEN << "inserted block(" << from << ", " << from << "):" << endl;
			cerr << *_blocks_pull[Association(from, from)] << endl << RESET;
		} else {
			cerr << BOLDYELLOW << "block(" << from << ", " << from << ") already exists:" << endl;
			cerr << *_blocks_pull[Association(from, from)] << endl << RESET;
		}

		f_it = _blocks_pull.find(Association(to, from));
		if(f_it == _blocks_pull.end()){
			_blocks_pull[Association(to, from)] = new SparseMatrixBlock();
			_blocks_pull[Association(to, from)]->setIdentity();
			cerr << BOLDGREEN << "inserted block(" << to << ", " << from << "):" << endl;
			cerr << *_blocks_pull[Association(to, from)] << endl << RESET;
		} else {
			cerr << BOLDYELLOW << "block(" << to << ", " << from << ") already exists:" << endl;
			cerr << *_blocks_pull[Association(to, from)] << endl << RESET;
		}

		//! PROBLEM HERE PORCO DIO:
		//! does not find block(54, 54) even if it exists -> creates new block but it already exists -> segfault
		f_it = _blocks_pull.find(Association(to, to));
		if(f_it == _blocks_pull.end()){
			_blocks_pull[Association(to, to)] = new SparseMatrixBlock();
			_blocks_pull[Association(to, to)]->setIdentity();
			cerr << BOLDGREEN << "inserted block(" << to << ", " << to << "):" << endl;
			cerr << *_blocks_pull[Association(to, to)] << endl << RESET;
		} else {
			cerr << BOLDYELLOW << "block(" << to << ", " << to << ") already exists:" << endl;
			cerr << *_blocks_pull[Association(to, to)] << endl << RESET;
		}
		cin.get();
	}

	_H = new SparseBlockMatrix(_vertices, _blocks_pull);

	//! Storage for U and L
	_L = _H->cholesky();
	_U = _L->transpose();

	//! Clean-up
	for (BlocksMap::iterator it = _blocks_pull.begin(); it != _blocks_pull.end(); ++it){
		it->second->setZero();
	}

	//! Storage for the rhs vector
	_B.init(_vertices.size());

	//! Storage for the increment vector
	_Y.init(_vertices.size());
	_dX.init(_vertices.size());
	return;
}

void SparseOptimizer::oneStep(void){
	real_ step_chi = 0.0;
	int step_inliers = 0;

	linearizeFactor(step_chi, step_inliers);
	cerr << BOLDWHITE << "inliers pose-pose = " << BOLDBLUE << step_inliers << "\t" << BOLDWHITE
			<< "chi pose-pose = " << BOLDBLUE<< step_chi << RESET << endl;

	//! Build and solve the linear system HdX = B;
	SparseMatrixBlock& h_00 = (*_blocks_pull[Association(0,0)]);
	h_00 += SparseMatrixBlock::Identity() * 1000000.0; //! Not good but ok for now

	_H->updateCholesky(_L);
	_L->updateTranspose(_U);

	DenseBlockVector y;
	_L->forwSub(_B, _Y);
	_U->backSub(_Y, _dX);

//	hessian->solveLinearSystem(_B, _dX);

	for (int i = 0; i < _dX.size; ++i) {
		Pose new_pose = Pose::Identity();
		new_pose = v2t(-(*_dX.blocks[i])) * _vertices[i].data();
		_vertices[i].setData(new_pose);
	}

	_dX.clear();
	_Y.clear();
	_B.clear();

	//! Clean-up the factors
	for (BlocksMap::iterator it = _blocks_pull.begin(); it != _blocks_pull.end(); ++it){
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
		SparseMatrixBlock& H_ii = (*_blocks_pull[Association(from, from)]);
		SparseMatrixBlock& H_ji = (*_blocks_pull[Association(to, from)]);
		SparseMatrixBlock& H_jj = (*_blocks_pull[Association(to, to)]);
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

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
  _jacobians_workspace.reset();

  _B.reset();
  _Y.reset();
  _dX.reset();
}

void SparseOptimizer::init(const VerticesContainer& vertices_,
    const EdgesContainer& edges_){
  _vertices = vertices_;
  _edges = edges_;

  //! TODO: solve problem here while creating storage -> segfault
  //! Initialize the Factors
  HessianBlocksMap::iterator f_it;
  for (int i = 0; i < _edges.size(); ++i) {
    VerticesContainer::const_iterator pose_i = std::find(_vertices.begin(),
        _vertices.end(), _edges[i].association().first);
    VerticesContainer::const_iterator pose_j = std::find(_vertices.begin(),
        _vertices.end(), _edges[i].association().second);
    int i_idx = pose_i->index();
    int j_idx = pose_j->index();

    _factors.push_back(Factor(i_idx, j_idx));
  }

  //! TODO
  //! FactorsVector ordered_factor = reorder(_factors, ordering_type);
  //! _jacobians_workspace.allocate(ordered_factors); -> also in the linearize function
  _jacobians_workspace.allocate(_factors);

  _H = SparseBlockMatrix(_vertices.size(), _jacobians_workspace);

  //! Allocate Matrices
  _H.allocateCholesky(_L);
  _L.allocateTransposed(_U);

  _jacobians_workspace.clear();

  //! Storage for the rhs vector
  _B.init(_vertices.size());

  //! Storage for the increment vector
  _Y.init(_vertices.size());
  _dX.init(_vertices.size());

  cerr << BOLDGREEN << "Init done" << endl << RESET;
  return;
}

void SparseOptimizer::oneStep(void){
  real_ step_chi = 0.0;
  int step_inliers = 0;

  linearizeFactor(step_chi, step_inliers);
  cerr << BOLDWHITE << "inliers pose-pose = " << BOLDBLUE << step_inliers << "\t" << BOLDWHITE
      << "chi pose-pose = " << BOLDBLUE<< step_chi << RESET << endl;

  //! Build and solve the linear system HdX = B;
  SparseMatrixBlock& h_00 = (*_jacobians_workspace.map.at(Association(0,0)));
  h_00 += SparseMatrixBlock::Identity() * 1000000.0; //! Not good but ok for now

  _H.computeCholesky(_L);
  _L.computeTranspose(_U);

  DenseBlockVector y;
  _L.forwSub(_B, _Y);
  _U.backSub(_Y, _dX);

  for (int i = 0; i < _dX.size; ++i) {
    Pose new_pose = Pose::Identity();
    new_pose = v2t(-(*_dX.blocks[i])) * _vertices[i].data();
    _vertices[i].setData(new_pose);
  }

  _dX.clear();
  _Y.clear();
  _B.clear();

  //! cleaning
  _jacobians_workspace.clear();
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

    int i_idx = pose_i->index();
    int j_idx = pose_j->index();

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
    SparseMatrixBlock& H_ii = (*_jacobians_workspace.map.at(Association(i_idx, i_idx)));
    SparseMatrixBlock& H_ji = (*_jacobians_workspace.map.at(Association(j_idx, i_idx)));
    SparseMatrixBlock& H_jj = (*_jacobians_workspace.map.at(Association(j_idx, j_idx)));

    H_ii.noalias() += Ji.transpose() * omega * Ji;
    H_ji.noalias() += Jj.transpose() * omega * Ji;
    H_jj.noalias() += Jj.transpose() * omega * Jj;

    //! Compute the factor contribution to the rhs vector
    DenseVectorBlock& b_i = (*_B.blocks[i_idx]);
    DenseVectorBlock& b_j = (*_B.blocks[j_idx]);
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

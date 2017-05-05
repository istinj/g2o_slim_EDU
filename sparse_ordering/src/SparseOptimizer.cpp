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
  _kernel_threshold = 100.0;
  _convergence_threshold = 1e-5;
  _num_iterations = 100;
  _total_chi = 0.0;
  _total_inliers = 0;
}

SparseOptimizer::~SparseOptimizer() {
  _jacobians_workspace.reset();

  _B.reset();
  _Y.reset();
  _dX.reset();
}

void SparseOptimizer::init(const VerticesVector& vertices_,
                           const EdgesContainer& edges_){
  cerr << BOLDWHITE << "Allocating memory for the optimizer" << endl << RESET;

  _vertices = vertices_;
  _edges = edges_;

  //! Create the Factors
  HessianBlocksMap::iterator f_it;
  for (int i = 0; i < _edges.size(); ++i) {
    VerticesVector::const_iterator pose_i = std::find(_vertices.begin(),
        _vertices.end(), _edges[i].association().first);
    VerticesVector::const_iterator pose_j = std::find(_vertices.begin(),
        _vertices.end(), _edges[i].association().second);
    int i_idx = pose_i->index();
    int j_idx = pose_j->index();

    _factors.push_back(Factor(i_idx, j_idx));
  }

  //! Allocate memory for the Hessian
  _jacobians_workspace.allocate(_factors);

  //! Create the sparseblockmatrix using the jacbian_workspace as memory
  bool exp = false;
  _H = SparseBlockMatrix(_vertices.size(), _jacobians_workspace);
  if (exp) {
    ofstream file_plain("../data/h_PLAIN.txt");
    file_plain << _H << endl;
    file_plain.close();
  }
  computeAMDPermutation(_hessian_permutation, _H);
  _H.reorder(_hessian_permutation, _jacobians_workspace, _factors);


  if (exp) {
    ofstream file_amd("../data/h_AMD.txt");
    file_amd << _H << endl;
    file_amd.close();
  }


  //! Allocate other matrices
  _H.allocateCholesky(_L);

  if (exp) {
    ofstream file_l("../data/chol_lower_AMD.txt");
    file_l << _L << endl;
    file_l.close();
  }

  _L.allocateTransposed(_U);

  if (exp) {
    ofstream file_u("../data/chol_upper_AMD.txt");
    file_u << _U << endl;
    file_u.close();
  }

  real_ fill_in_H = (_vertices.size()+(_edges.size()*2))/(real_)(_vertices.size()*_vertices.size());
  real_ rate = _L.nnz()/(real_)(_H.nnz() - _factors.size());
  cerr << BOLDWHITE << "Hessian fill-in rate: \t" << BOLDGREEN << fill_in_H << endl;
  cerr << BOLDWHITE << "Cholesky fill-in rate: \t" << BOLDGREEN << rate * fill_in_H << endl << endl;

  _jacobians_workspace.clear();

  //! Storage for the rhs vector
  _B.init(_vertices.size());

  //! Storage for the increment vector
  _Y.init(_vertices.size());
  _dX.init(_vertices.size());
  return;
}


void SparseOptimizer::converge(void) {
  real_ prev_total_chi = 0.0;
  int iter_cnt = 0;
  bool suppress_outliers = false;

  for (int i = 0; i < _num_iterations; ++i) {
    oneStep(suppress_outliers);
    ++iter_cnt;

    const real_ delta_total_chi = fabs(prev_total_chi - _total_chi);
    if (_kernel_threshold > delta_total_chi) {
      suppress_outliers = true;
      oneStep(suppress_outliers);
      oneStep(suppress_outliers);
      oneStep(suppress_outliers);
      iter_cnt += 3;
      break;
    } else {
      prev_total_chi = _total_chi;
    }

    if ( i == _num_iterations - 1) {
      cerr << BOLDRED << "System did not converged| total error: " << delta_total_chi << RESET << endl;
      return;
    }
  }
  cerr << BOLDGREEN << "System converged in " << iter_cnt << " iterations" << RESET << endl;
}

void SparseOptimizer::oneStep(bool suppress_outliers_){
  std::chrono::high_resolution_clock::time_point t_0, t_1;
  t_0 = std::chrono::high_resolution_clock::now();
  _total_chi = 0.0;
  _total_inliers = 0;

  linearizeFactor(_total_chi, _total_inliers, suppress_outliers_);
  cerr << BOLDWHITE << "inliers pose-pose = " << BOLDBLUE << _total_inliers << "\t" << BOLDWHITE
      << "chi pose-pose = " << BOLDBLUE<< _total_chi << RESET << endl;

  //! Build and solve the linear system HdX = B;
  _jacobians_workspace(0,0) += SparseMatrixBlock::Identity() * 1000000.0; //! Not good but ok for now

  _U.printBlock(2,2);
  _H.computeCholesky(_L);
  _L.computeTranspose(_U);

  DenseBlockVector y;
  _L.forwSub(_B, _Y);
  _U.backSub(_Y, _dX);

  updateVertices();

  //! cleaning
  _dX.clear();
  _Y.clear();
  _B.clear();

  _jacobians_workspace.clear();
  t_1 = std::chrono::high_resolution_clock::now();
  double execution_time = (std::chrono::duration_cast<std::chrono::microseconds>(t_1 - t_0).count() / 1e06);
  std::cerr << BOLDWHITE << "Step time:\t" << BOLDGREEN << execution_time << "s" << std::endl << RESET;
}


void SparseOptimizer::updateGraph(Graph& graph_) {
  if (graph_.vertices().size() != _vertices.size())
    throw std::runtime_error("Error, number of vertices must remain the same");
  graph_.updateVertices(_vertices);
}

void SparseOptimizer::computeAMDPermutation(IntVector& permutation_AMD_,
                                            SparseBlockMatrix& matrix_) {
  cs* cs_matrix = matrix_.toCs();
  int* amd_odering = cs_amd(1,cs_matrix);

  permutation_AMD_.resize(matrix_.numRows());

  for (int i = 0; i < matrix_.numRows(); ++i) {
    permutation_AMD_[amd_odering[i]] = i;
  }

  cs_spfree(cs_matrix);
  cs_free(amd_odering);
}


void SparseOptimizer::updateVertices(void) {
  Pose new_pose;
  for (int i = 0; i < _dX.size; ++i) {
    new_pose.setIdentity();
    new_pose = v2t(-(*_dX.blocks[i])) * _vertices[i].data();
    _vertices[i].setData(new_pose);
  }
}


void SparseOptimizer::linearizeFactor(real_& total_chi_,
                                      int& inliers_,
                                      bool suppress_outliers_){
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
    VerticesVector::const_iterator pose_i = std::find(_vertices.begin(),
        _vertices.end(), edge.association().first);
    VerticesVector::const_iterator pose_j = std::find(_vertices.begin(),
        _vertices.end(), edge.association().second);

    int i_idx = pose_i->index();
    int j_idx = pose_j->index();

    //! Compute error and jacobian
    errorAndJacobian(pose_i->data(), pose_j->data(), _edges[i].data(),
        e, Ji, Jj);

    //! Robust kernel
    real_ chi = e.transpose() * omega * e;
    if(chi > _kernel_threshold) {
      if (suppress_outliers_)
        continue;
      omega *= sqrtf(_kernel_threshold / chi);
      chi = _kernel_threshold;
    } else {
      inliers_++;
    }
    total_chi_ += chi;

    //! Compute the factor contribution to the Hessian
    SparseMatrixBlock& H_ii = _jacobians_workspace(i_idx, i_idx);
    SparseMatrixBlock& H_ji = _jacobians_workspace(j_idx, i_idx);
    SparseMatrixBlock& H_jj = _jacobians_workspace(j_idx, j_idx);

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


void SparseOptimizer::errorAndJacobian(const Pose& xi,
                                       const Pose& xj,
                                       const PoseMeas& zr,
                                       Vector12& error,
                                       Matrix12_6& Ji,
                                       Matrix12_6& Jj) {
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

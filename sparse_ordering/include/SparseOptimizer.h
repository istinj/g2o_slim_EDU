/*
 * SparseOptimizer.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef SPARSEOPTIMIZER_H_
#define SPARSEOPTIMIZER_H_

#include "defs.h"
#include "Graph.h"
#include "Vertex.h"
#include "Edge.h"
#include "Factor.h"
#include "Workspace.h"
#include "SparseBlockMatrix.h"

namespace sparse {

typedef std::vector<sparse::Vertex> VerticesVector;
typedef std::vector<sparse::Edge> EdgesContainer;
typedef std::vector<Factor> FactorsVector;

class SparseOptimizer {
public:
  SparseOptimizer();
  virtual ~SparseOptimizer();

  //! This will allocate all the memory that will contain the Hessian, its Cholesky and the transposed
  //! of the Cholesky.
  void init(const VerticesVector& vertices_,
            const EdgesContainer& edges_,
            const SolverType& type_);
  void converge(void);
  void updateGraph(Graph& graph_);

protected:
  void _computeAMDPermutation(IntVector& permutation_AMD_,
                             SparseBlockMatrix& matrix_);
  //! This function acutally computes all the stuff allocated during the INIT function
  void _oneStep(bool suppress_outliers_);
  void _updateVertices(void);
  void _linearizeFactor(Real& total_chi,
                       int& inliers_,
                       bool suppress_outliers_);
  void _errorAndJacobian(const Pose& xi,
                        const Pose& xj,
                        const PoseMeas& zr,
                        Vector12& error,
                        Matrix12_6& Ji,
                        Matrix12_6& Jj);

  VerticesVector _vertices;
  EdgesContainer _edges;

  Workspace _jacobians_workspace;
  FactorsVector _factors;
  IntVector _permutation;

  DenseBlockVector _B;
  DenseBlockVector _Y;
  DenseBlockVector _dX;

  SparseBlockMatrix _H;
  SparseBlockMatrix _L;
  SparseBlockMatrix _U;

  Real _kernel_threshold;
  Real _convergence_threshold;
  Real _total_chi;
  int _total_inliers;
  int _num_iterations;

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
};

} /* namespace sparse */

#endif /* SPARSEOPTIMIZER_H_ */

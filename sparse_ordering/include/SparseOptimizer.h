/*
 * SparseOptimizer.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef SPARSEOPTIMIZER_H_
#define SPARSEOPTIMIZER_H_

#include "defs.h"
#include "Vertex.h"
#include "Edge.h"
#include "Factor.h"
#include "SparseBlockMatrix.h"

namespace sparse {

typedef std::vector<sparse::Vertex> VerticesContainer;
typedef std::vector<sparse::Edge> EdgesContainer;
typedef std::map<Factor, SparseMatrixBlock*, FactorComparator> FactorsMap;

class SparseOptimizer {
public:
	SparseOptimizer();
	virtual ~SparseOptimizer();

	void init(const VerticesContainer& vertices_,
			const EdgesContainer& edges_);
	void oneStep(void);

protected:
	void linearizeFactor(real_& total_chi, int& inliers_);
	void errorAndJacobian(const Pose& xi, const Pose& xj,const PoseMeas& zr,
			Vector12& error, Matrix12_6& Ji, Matrix12_6& Jj);

	VerticesContainer _vertices;
	EdgesContainer _edges;

	//! TODO: 	is this container for the matrices good? Do i need something more complex?
	FactorsMap _factors; //! Matrices Pull -> not factors, name must be refactored

	DenseBlockVector _B;
	DenseBlockVector _Y;
	DenseBlockVector _dX;

	SparseBlockMatrix* _H;
	SparseBlockMatrix* _L;
	SparseBlockMatrix* _U;

	real_ _kernel_threshold = 1000.0;
	//! TODO:	How to generate orderigs? How is structured the vector of ints containing
	//! 		the ordering??

	//! TODO	Initialize memory for L and U in the init function.

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

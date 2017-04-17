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

namespace sparse {

struct Factor {
	int from;
	int to;

	SparseMatrixBlock* block_ii;
	SparseMatrixBlock* block_ji;
	SparseMatrixBlock* block_jj;

	inline void clear(void) {
		from = 0;
		to = 0;

		block_ii->setZero();
		block_ji->setZero();
		block_jj->setZero();
	}

	//! TODO: use shared ptr
//	std::shared_ptr<SparseMatrixBlock> block_11;
//	std::shared_ptr<SparseMatrixBlock> block_21;
//	std::shared_ptr<SparseMatrixBlock> block_22;
};

struct DenseBlockVector {
	int num_block;
	DenseVectorContainer blocks;

	inline void clear(void) {
		for (int i = 0; i < blocks.size(); ++i) {
			blocks[i]->setZero();
		}
	}

	inline void init(const int size_){
		blocks.resize(size_);
		num_block = size_;
	}
};

typedef std::vector<sparse::Vertex> VerticesContainer;
typedef std::vector<sparse::Edge> EdgesContainer;
typedef std::vector<sparse::Factor> FactorsContainer;

class SparseOptimizer {
public:
	SparseOptimizer();
	virtual ~SparseOptimizer();

	void init(const VerticesContainer& vertices_,
			const EdgesContainer& edges_);
	void oneStep(void);

protected:
	//! TODO: container + SparseBlockMatrix + DenseVector + LinearSolver
	void linearizeFactor(real_& total_chi, int& inliers_);
	void errorAndJacobian(const Pose& xi, const Pose& xj,const PoseMeas& zr,
			Vector12& error, Matrix12_6& Ji, Matrix12_6& Jj);

	VerticesContainer _vertices;
	EdgesContainer _edges;
	FactorsContainer _factors;

	DenseBlockVector _B;

	real_ _kernel_threshold = 1000.0;

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

/*
 * Workspace.h
 *
 *  Created on: 28/apr/2017
 *      Author: istin
 */

#ifndef WORKSPACE_H_
#define WORKSPACE_H_

#include <defs.h>
#include <factor.h>

namespace sparse {

typedef IndicesSparseMatrixBlocksPtrMap WorkspaceMap;

struct Workspace {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Workspace();
  virtual ~Workspace();

  //! DESTROYS the memory allocated
  void reset(void);
  //! Sets to 0 all the blocks allocated
  void clear(void);
  //! ALLOCATES blocks and sets them to identity
  void allocate(const std::vector<Factor>& factors_);
  void allocate(const std::vector<IntPair>& indices_);
  void allocateOneBlock(const IntPair& indices_);
  void allocateOneBlock(const int r_, const int c_);

  SparseMatrixBlock& operator()(const IntPair& indices_) const;
  SparseMatrixBlock& operator()(const int r_, const int c_) const;

  //! TODO c-style Matrix of SparseMatrixBlock* -> SparseMatrixBlock***
  WorkspaceMap memory_map;
  int dimension;

};

} /* namespace sparse */

#endif /* WORKSPACE_H_ */

/*
 * Workspace.h
 *
 *  Created on: 28/apr/2017
 *      Author: istin
 */

#ifndef WORKSPACE_H_
#define WORKSPACE_H_

#include <defs.h>
#include <Factor.h>

namespace sparse {

typedef BlocksMap WorkspaceMap;
typedef BlocksMultiMap WorkspaceMultiMap;

struct Workspace {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Workspace();
  virtual ~Workspace();

  void reset(void);
  void clear(void);
  void allocate(std::vector<Factor>& factors_);
  void allocate(std::vector<Association>& indices_);

  //! TODO c-style Matrix of SparseMatrixBlock*
  WorkspaceMap map;
  int dimension;

};

} /* namespace sparse */

#endif /* WORKSPACE_H_ */

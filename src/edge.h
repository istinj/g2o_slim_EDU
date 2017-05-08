/*
 * Edge.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef EDGE_H_
#define EDGE_H_

#include "defs.h"

namespace sparse {

class Edge {
public:
  Edge();
  virtual ~Edge();

  void setEdge(const IntPair& association_,
      const PoseMeas& data_,
      const OmegaPose& omega_);
  void setEdge(const IntPair& association_,
      const int& sensor_ID_,
      const PoseMeas& data_,
      const OmegaPose& omega_);

  inline const int sensorID(void) const {return _sensor_id;};
  inline const IntPair& association(void) const {return _id_association;};
  inline const PoseMeas& data(void) const {return _data;};
  inline const OmegaPose& omega(void) const {return _omega;};

protected:
  IntPair _id_association;
  int _sensor_id;

  PoseMeas _data;
  OmegaPose _omega;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* EDGE_H_ */

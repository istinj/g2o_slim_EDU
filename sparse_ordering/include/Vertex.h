/*
 * Vertex.h
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#ifndef VERTEX_H_
#define VERTEX_H_

#include "defs.h"

namespace sparse {

class Vertex {
public:
	Vertex();
	Vertex(const int& id_, const int& index_, const Pose& data_);
	virtual ~Vertex();

	void setVertex(const int& id_, const int& index_, const Pose& data_);
	void setData(const Pose& data_);

	inline bool operator==(const Vertex& other_) const {return this->_id == other_._id;};
	inline bool operator==(const int& ID) const {return this->_id == ID;};

	inline const Pose& data(void) const {return _data;};
	inline const int id(void) const {return _id;};
	inline const int index(void) const {return _index;};

protected:
	Pose _data;
	int _id;
	int _index;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

} /* namespace sparse */

#endif /* VERTEX_H_ */

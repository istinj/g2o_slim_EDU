/*
 * Vertex.h
 *
 *  Created on: 06/mar/2017
 *      Author: istin
 */

#ifndef VERTEX_H_
#define VERTEX_H_

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace optimizer {

typedef Eigen::Isometry3f Pose;
typedef Eigen::Vector3f LandmarkXYZ;

template <class _DataType>
class Vertex{
public:
	Vertex();
	Vertex(const int id_, const _DataType& data_);
	Vertex(const int id_, const _DataType& data_, const int index_);
	~Vertex();

	inline bool operator==(const Vertex<_DataType>& other_) const {return this->_id == other_._id;};
	inline bool operator==(const int ID) const {return this->_id == ID;};

	void setVertex(const int id_, const _DataType& data_);
	void setVertex(const int id_, const _DataType& data_, const int index_);
	void setData(const _DataType& data_);

	inline const _DataType& data(void) const {return _data;};
	inline const int id(void) const {return _id;};
	inline const int index(void) const {return _index;};

private:
	_DataType _data;
	int _id;
	int _index = -1; //index of the vertices in its container

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <class _DataType>
Vertex<_DataType>::Vertex(){
	_id = -1;
}

template <class _DataType>
Vertex<_DataType>::Vertex(const int id_, const _DataType& data_){
	_id = id_;
	_data = data_;
}

template <class _DataType>
Vertex<_DataType>::Vertex(const int id_, const _DataType& data_, const int index_){
	_id = id_;
	_data = data_;
	_index = index_;
}


template <class _DataType>
Vertex<_DataType>::~Vertex(){
}


template <class _DataType>
void Vertex<_DataType>::setVertex(const int id_, const _DataType& data_){
	_id = id_;
	_data = data_;
}

template <class _DataType>
void Vertex<_DataType>::setVertex(const int id_, const _DataType& data_, const int index_){
	_id = id_;
	_data = data_;
	_index = index_;
}

template <class _DataType>
void Vertex<_DataType>::setData(const _DataType& data_){
   _data = data_;
}

}/* namespace optimizer */

#endif /* VERTEX_H_ */

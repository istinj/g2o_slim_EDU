/*
 * Edge.h
 *
 *  Created on: 06/mar/2017
 *      Author: istin
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace optimizer {

typedef Eigen::Vector3f PointMeas;
typedef Eigen::Matrix3f OmegaPoint;
typedef Eigen::Isometry3f PoseMeas;
typedef Eigen::Matrix<float, 6, 6> OmegaPose;

template <class _DataType, class _InfoMatrixType>
class Edge {
public:
	Edge();
	~Edge();

	void setEdge(const std::pair<int, int> association_,
			const _DataType data_,
			const _InfoMatrixType omega_);

	void setEdge(const std::pair<int, int> association_,
			const int sensor_ID_,
			const _DataType data_,
			const _InfoMatrixType omega_);

	inline const int sensorID(void) const {return _sensor_id;};
	inline const std::pair<int, int> association(void) const {return _IDassociation;};
	inline const _DataType& data(void) const {return _data;};
	inline const _InfoMatrixType& omega(void) const {return _Omega;};

private:
	std::pair<int, int> _IDassociation;
	int _sensor_id;

	_DataType _data;
	_InfoMatrixType _Omega;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <class _DataType, class _InfoMatrixType>
Edge<_DataType, _InfoMatrixType>::Edge(){
	_IDassociation = std::pair<int, int>(-1, -1);
	_sensor_id = -1; //! also for odometry meas
}

template <class _DataType, class _InfoMatrixType>
Edge<_DataType, _InfoMatrixType>::~Edge(){
}

template <class _DataType, class _InfoMatrixType>
void Edge<_DataType, _InfoMatrixType>::setEdge(const std::pair<int, int> association_,
		const _DataType data_,
		const _InfoMatrixType omega_){
	_IDassociation = association_;
	_data = data_;
	_Omega = omega_;
}

template <class _DataType, class _InfoMatrixType>
void Edge<_DataType, _InfoMatrixType>::setEdge(const std::pair<int, int> association_,
		const int sensor_ID_,
		const _DataType data_,
		const _InfoMatrixType omega_){
	_IDassociation = association_;
	_sensor_id = sensor_ID_;
	_data = data_;
	_Omega = omega_;
}
}/* namespace optimizer */

#endif /* EDGE_H_ */

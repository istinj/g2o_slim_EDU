/*
 * Edge.cpp
 *
 *  Created on: 11/apr/2017
 *      Author: istin
 */

#include <Edge.h>

namespace sparse {

Edge::Edge() {
	// TODO Auto-generated constructor stub
	_id_association = Association(-1, -1);
	_sensor_id = -1;
}

Edge::~Edge() {
	// TODO Auto-generated destructor stub
}

void Edge::setEdge(const Association& association_,
		const PoseMeas& data_,
		const OmegaPose& omega_){
	_id_association = association_;
	_data = data_;
	_omega = omega_;
}
void Edge::setEdge(const Association& association_,
		const int& sensor_ID_,
		const PoseMeas& data_,
		const OmegaPose& omega_){
	_id_association = association_;
	_sensor_id = sensor_ID_;
	_data = data_;
	_omega = omega_;
}

} /* namespace sparse */

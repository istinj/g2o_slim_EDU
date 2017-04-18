/*
 * Factor.cpp
 *
 *  Created on: 17/apr/2017
 *      Author: istin
 */

#include <Factor.h>

namespace sparse {
void Factor::clear(void) {
	from = 0;
	to = 0;

//	block_ii->setZero();
//	block_ji->setZero();
//	block_jj->setZero();
}

bool Factor::operator ==(const Factor& other_) const {
	if (this->from == other_.from && this->to == other_.to)
		return true;
	return false;
}

bool Factor::operator ==(const Association& factor_association_) const {
	if (this->from == factor_association_.first && this->to == factor_association_.second)
		return true;
	return false;
}

} /* namespace sparse */

/*
 * NBoxFBC.cpp
 *
 *  Created on: 13.07.2017
 *      Author: ciesla
 */

#include "NBoxFBC.h"

NBoxFBC::NBoxFBC(double s, double ndx, [[maybe_unused]] double vdx) : Surface(s, ndx) {
	// TODO Auto-generated constructor stub
}

double NBoxFBC::getArea() const {
	return pow(this->size, 2);
}

RSAVector NBoxFBC::getTranslation([[maybe_unused]] double s, [[maybe_unused]] const RSAVector &p1,
                                  [[maybe_unused]] const RSAVector &p2)
{
	return {};
}

RSAVector NBoxFBC::getTranslation(const RSAVector &p1, const RSAVector &p2) const {
	return NBoxFBC::getTranslation(this->size, p1, p2); }

RSAVector NBoxFBC::vector(const RSAVector &v) const {
	return this->vectorFreeBC(v);
}

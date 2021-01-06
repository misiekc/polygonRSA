/*
 * BoundaryConditions.h
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_

#include "Vector.h"

class BoundaryConditions {
public:
    using RSAVector = Vector<2>;

	virtual ~BoundaryConditions() = default;

	virtual double distance2(const RSAVector &p1, const RSAVector &p2) const = 0;

	/**
	 * @brief Returns translation that should be applied to @a p2 to move him to the "proximity" of @a p1
	 */
	virtual RSAVector getTranslation(const RSAVector& p1, const RSAVector &p2) const = 0;

};

#endif /* BOUNDARYCONDITIONS_H_ */

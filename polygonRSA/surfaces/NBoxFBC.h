/*
 * NBoxFBC.h
 *
 *  Created on: 13.07.2017
 *      Author: ciesla
 */

#ifndef SURFACES_NBOXFBC_H_
#define SURFACES_NBOXFBC_H_

#include "../Surface.h"
#include "../RND.h"

class NBoxFBC: public Surface {
private:
	static RSAVector getTranslation(double s, const RSAVector &p1, const RSAVector &p2);

public:
	NBoxFBC(double s, double ndx, double vdx);

	double getArea() const override;
	RSAVector getTranslation(const RSAVector &p1, const RSAVector &p2) const override;
	RSAVector vector(const RSAVector &v) const override;
};

#endif /* SURFACES_NBOXFBC_H_ */

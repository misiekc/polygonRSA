/*
 * NBoxPBC.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACES_NBOXPBC_H_
#define SURFACES_NBOXPBC_H_

#include "../Surface.h"
#include "../RND.h"

class NBoxPBC : public Surface {
private:
	static RSAVector getTranslation(double s, const RSAVector &p1, const RSAVector &p2);

public:
	NBoxPBC(double s, double ndx, double vdx);

	double getArea() const override;
	RSAVector getTranslation(const RSAVector &p1, const RSAVector &p2) const override;
	RSAVector vector(const RSAVector &v) const override;

	RSAVector checkPosition(const RSAVector &da) const override;
};

#endif /* SURFACES_NBOXPBC_H_ */

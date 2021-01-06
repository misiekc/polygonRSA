/*
 * Surface.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef SURFACE_H_
#define SURFACE_H_


#include <vector>
//#include "DeviceNeighbourGrid.cu"

#include "shape/Shape.h"
#include "NeighbourGrid.h"

class Surface : public BoundaryConditions {

protected:

	NeighbourGrid<const Shape> *list;

	double size;

	RSAVector vectorFreeBC(const RSAVector &v) const;
	RSAVector vectorPeriodicBC(const RSAVector &v) const;

public:
	Surface(double s, double ndx);
	virtual ~Surface();

	void clear();
	void add(const Shape *s);
	const Shape *check(const Shape *s);
	void getNeighbours(std::vector<const Shape*> *result, const RSAVector &da);
	const Shape *getClosestNeighbour(const RSAVector &da);
	const Shape *getClosestNeighbour(const RSAVector &da, const std::vector<const Shape*> &neighbours);
	NeighbourGrid<const Shape> *getNeighbourGrid();
    double distance2(const RSAVector &a1, const RSAVector &a2) const override;

	RSAVector getTranslation(const RSAVector &p1, const RSAVector &p2) const override = 0;
	virtual RSAVector vector(const RSAVector &v) const = 0;
	virtual RSAVector checkPosition(const RSAVector &da) const;
	virtual double getArea() const = 0;

//	void drawShapes(Graphics g, double scale);
//	void drawShapes(Graphics g, double scale, double[] ta);
};

#endif /* SURFACE_H_ */

/*
 * RoundedPolygon.h
 *
 *  Created on: 27.03.2020
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_ROUNDEDPOLYGON_H_
#define SHAPES_POLYGONS_ROUNDEDPOLYGON_H_

#include <vector>
#include <utility>
#include <cstddef>

#include "../../AnisotropicShape2D.h"
#include "../../../Vector.h"

#include "Polygon.h"

class RoundedPolygonAlgebraic : public Polygon {

private:
	static double distance2pq(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q);
	static void gradientDistance2pq(std::array<double, 2> &gradient, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q);
	static void normalizeVolume(std::istringstream &in);
	static double getArea();


protected:
	static double radius;

    static double calculateCircumscribedCircleRadius();

    static double lineLineDistance2(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
    static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

	bool overlapComplexCheck(RSAVector &position, double angle, RSAVector &polposition, double polangle) const override;
    bool voxelInsideComplexCheck(const RSAVector &spatialCenter, double halfSpatialSize, double angularCenter, double halfAngularSize) const override;

public:
//    RoundedPolygon();

	static void initClass(const std::string &args);

	Shape *clone() const override;

//	bool overlap(BoundaryConditions<2> *bc, const RSAShape *s) const override;

//	bool voxelInside(BoundaryConditions<2> *bc, const RSAVector &voxelPosition, const Orientation<1> &voxelOrientation,

	double getVolume() const override;

	std::string toPovray() const override;
	std::string toWolfram() const override;
};

#endif /* SHAPES_POLYGONS_ROUNDEDPOLYGON_H_ */

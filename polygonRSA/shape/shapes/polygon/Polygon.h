/*
 * Polygon.h
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#ifndef SHAPES_POLYGONS_POLYGON_H_
#define SHAPES_POLYGONS_POLYGON_H_

#include <vector>
#include <utility>
#include <cstddef>

#include "../../AnisotropicShape2D.h"
#include "../../../Vector.h"

class Polygon : public RSAShape {
private:
    static void normalizeVolume(std::istringstream &in);
    static bool pointInsidePolygon(const RSAVector &point, const std::vector<RSAVector> &vertiecs);

    bool voxelInsideFullAngleCheck(const RSAVector &spatialCenter, double halfSpatialSize) const;
    bool pointInsidePushedVertices(const RSAVector &point, double pushDistance) const;
    bool pointInsidePushedEdges(const RSAVector &point, double pushDistance) const;

    static constexpr std::size_t INSPHERE_SEARCH_DIVISIONS = 20;
    static constexpr double INSPHERE_SEARCH_FACTOR = M_SQRT2;
    static constexpr double INSPHERE_SEARCH_PRECISION = 1e-8;

protected:
    //polar coordinates of all vertices
    static std::vector<double> vertexR;
    static std::vector<double> vertexTheta;
    static std::vector<std::pair<size_t, size_t>> segments;
    static std::vector<std::pair<size_t, size_t>> helperSegments;

    static void clearOldData();
    static void parseVertices(std::istringstream &in);
    static void parseSegments(std::istringstream &in);
    static void parseHelperSegments(std::istringstream &in);

    static RSAVector getStaticVertexPosition(std::size_t index);
    static std::vector<RSAVector> getStaticVerticesPositions();

    static double calculateCircumscribedCircleRadius();
    static double calculateInscribedCircleRadius(const RSAVector &origin = RSAVector{});
    static void centerPolygon();
    static void createStarHelperSegments();

    //calculate the area of the triangle made from the origin, vertex i, and vertex j
    static double getTriangleArea(size_t i, size_t j);

    //test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
    static bool lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

    //same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
    static bool lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4);

    static double segmentPointDistance2(const RSAVector &s1, const RSAVector &s2, const RSAVector &point);

    virtual bool overlapComplexCheck(RSAVector &position, double angle, RSAVector &polposition, double polangle) const;
    virtual bool voxelInsideComplexCheck(const RSAVector &spatialCenter, double halfSpatialSize, double angularCenter,
                                         double halfAngularSize) const;

    RSAVector getVertexPosition(std::size_t index) const;
    std::vector<RSAVector> getVerticesPositions() const;

    void vertexToPovray(std::size_t index, std::ostream &out) const;

#ifdef CUDA_ENABLED
    static void cuInit();
	static void cuFree();
#endif


public:
    static void initClass(const std::string &args);

    static const std::vector<double> &getVertexR() { return vertexR; }
    static const std::vector<double> &getVertexTheta() { return vertexTheta; }
    static const std::vector<std::pair<size_t, size_t>> &getSegments() { return segments; }
    static const std::vector<std::pair<size_t, size_t>> &getHelperSegments() { return helperSegments; }

    RSAShape *clone() const override;
    double getVolume() const override;

    bool overlap(RSABoundaryConditions *bc, const RSAShape *s) const override;

    bool voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition, const RSAOrientation &voxelOrientation,
                     double spatialSize, double angularSize) const override;
    std::string toPovray() const override;
    std::string toString() const override;
    std::string toWolfram() const override;

    static bool isPolygonConvex();

    static RSAVector calculateOptimalOrigin();
};

#endif /* SHAPES_POLYGONS_POLYGON_H_ */
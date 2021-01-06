//
// Created by pkua on 06.04.2020.
//

#include "RegularRoundedPolygonAlgebraic.h"
#include "RegularRoundedPolygonGeometric.h"

void RegularRoundedPolygonAlgebraic::initClass(const std::string &attr) {
    RegularDiskopolygonAttributes attributes(attr);

    double nSides = attributes.getNSides();
    double radius = attributes.getRadius();
    double halfDiagonal = attributes.getHalfDiagonal();

    std::ostringstream roundedPolygonAttributes;
    roundedPolygonAttributes.precision(std::numeric_limits< double >::max_digits10);

    roundedPolygonAttributes << radius << " " << nSides << " rt ";

    for (std::size_t i{}; i < nSides; i++)
        roundedPolygonAttributes << " " << halfDiagonal << " " << (2*M_PI*i/nSides);

    roundedPolygonAttributes << " " << nSides;
    for (std::size_t i{}; i < nSides; i++)
        roundedPolygonAttributes << " " << i;

    RoundedPolygonAlgebraic::initClass(roundedPolygonAttributes.str());

    // Make angular voxel size smaller - regular polygon has an n-fold rotational symmetry
    ShapeStaticInfo shapeInfo = Shape::getShapeStaticInfo();
    shapeInfo.setAngularVoxelSize(2*M_PI/nSides);
    Shape::setShapeStaticInfo(shapeInfo);
}

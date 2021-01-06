//
// Created by Piotr Kubala on 26/03/2020.
//

#include <sstream>
#include <cmath>

#include "RegularDiskopolygon.h"
#include "..//SpheroCylinder2D.h"

namespace {
    class FreeBC : public RSABoundaryConditions {
    public:
        double distance2(const RSAVector &p1, const RSAVector &p2) const { return (p1 - p2).norm2(); }
        RSAVector getTranslation([[maybe_unused]]const RSAVector &p1, [[maybe_unused]] const RSAVector &p2) const {
            return RSAVector{};
        }
    };
}

RegularDiskopolygonAttributes RegularDiskopolygon::attributes;

void RegularDiskopolygon::initClass(const std::string &attr) {
    attributes = RegularDiskopolygonAttributes(attr);

    std::ostringstream spherocylinderAttr;
    spherocylinderAttr << getSideLength() << " " << getRadius();
    SpheroCylinder2D::initClass(spherocylinderAttr.str());

    ShapeStaticInfo shapeInfo;
    shapeInfo.setCircumsphereRadius(getHalfDiagonal() + getRadius());
    shapeInfo.setInsphereRadius(getHeight() + getRadius());
    shapeInfo.setAngularVoxelSize(2*M_PI / getNSides());
    shapeInfo.setSupportsSaturation(true);
    shapeInfo.setDefaultCreateShapeImpl<RegularDiskopolygon>();
    RSAShape::setShapeStaticInfo(shapeInfo);
}

void RegularDiskopolygonAttributes::normalizeVolume() {
    double volume = 0.25 * this->nSides * std::pow(this->sideLength, 2) / tan(M_PI / this->nSides)   // Polygon inside
                    + this->nSides * this->radius * this->sideLength                                  // Pushed edges
                    + M_PI * std::pow(this->radius, 2);                                   // Corners

    this->sideLength /= std::sqrt(volume);
    this->radius /= std::sqrt(volume);
    this->height /= std::sqrt(volume);
    this->halfDiagonal /= std::sqrt(volume);
}

bool RegularDiskopolygon::overlap(RSABoundaryConditions *bc, const RSAShape *s) const {
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    RegularDiskopolygon other = dynamic_cast<const RegularDiskopolygon &>(*s);
    this->applyBC(bc, &other);

    FreeBC freeBC;
    for (std::size_t i{}; i < getNSides(); i++) {
        SpheroCylinder2D sc1 = this->getSpherocylinder(i);
        for (std::size_t j{}; j < getNSides(); j++) {
            SpheroCylinder2D sc2 = other.getSpherocylinder(j);
            if (sc1.overlap(&freeBC, &sc2))
                return true;
        }
    }

    return false;
}

bool RegularDiskopolygon::voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition,
                                      const RSAOrientation &orientation, double spatialSize, double angularSize) const
{
    switch (this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    RSAVector defaultSpherocylinderOffset{{0, getHeight()}};

    for (std::size_t i{}; i < getNSides(); i++) {
        double virtualScAngleFrom = orientation[0] + static_cast<double>(i) * 2 * M_PI / getNSides();
        double virtualScAngleTo = virtualScAngleFrom + angularSize;
        AnisotropicShape2D::normalizeAngleRange(0, &virtualScAngleFrom, &virtualScAngleTo, 2*M_PI);

        RectangularBounding rectangularBounding = RectangularBoundingBuilder::forArch(defaultSpherocylinderOffset,
                                                                                      virtualScAngleFrom,
                                                                                      virtualScAngleTo);

        RSAVector voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
        rectangularBounding.expand(spatialSize);
        rectangularBounding.translate(voxelPosition + voxelTranslation);

        for (std::size_t j{}; j < getNSides(); j++) {
            FreeBC freeBC;

            auto sc = this->getSpherocylinder(j);
            if (sc.pointInside(&freeBC, rectangularBounding.getTopLeft(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getTopRight(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getBottomLeft(), virtualScAngleFrom, virtualScAngleTo) &&
                sc.pointInside(&freeBC, rectangularBounding.getBottomRight(), virtualScAngleFrom, virtualScAngleTo))
            {
                return true;
            }
        }
    }

    return false;
}

RSAShape *RegularDiskopolygon::clone() const {
    return new RegularDiskopolygon(*this);
}

SpheroCylinder2D RegularDiskopolygon::getSpherocylinder(std::size_t index) const {
    SpheroCylinder2D sc;

    double scAngle = this->getAngle() + static_cast<double>(index) * 2 * M_PI / getNSides();
    RSAVector defaultSpherocylinderOffset{{0, getHeight()}};
    sc.rotate({{scAngle}});
    sc.translate(this->getPosition() + Matrix<2, 2>::rotation(scAngle) * defaultSpherocylinderOffset);

    return sc;
}

std::string RegularDiskopolygon::toWolfram() const {
    std::ostringstream out;

    out << "{";
    for (std::size_t i{}; i < getNSides(); i++)
        out << this->getSpherocylinder(i).toWolfram() << ", ";
    // For this->getAngle() == 0, we want top side to be horizontal
    double mathematicaAngle = this->getAngle() + M_PI/2 + M_PI/getNSides();
    out << "RegularPolygon[" << this->getPosition() << ", {" << getHalfDiagonal() << ", " << mathematicaAngle << "}, ";
    out << getNSides() << "]}";

    return out.str();
}

void RegularDiskopolygon::setOrientation(const RSAOrientation &orientation) {
    double angle = AnisotropicShape2D::normalizeAngle(orientation[0], 2*M_PI/getNSides());
    RSAShape::setOrientation({{angle}});
}

void RectangularBoundingBuilder::recurseQuarters(RectangularBounding &bounding, const RSAVector &zeroAngleVector,
                                                 double angleTo, double quarterAngle)
{
    if (angleTo < quarterAngle) {
        bounding.addPoint(Matrix<2, 2>::rotation(angleTo) * zeroAngleVector);
    } else {
        bounding.addPoint(Matrix<2, 2>::rotation(quarterAngle) * zeroAngleVector);
        recurseQuarters(bounding, zeroAngleVector, angleTo, quarterAngle + M_PI / 2);
    }
}

void RectangularBounding::addPoint(const RSAVector &p) {
    if (p[0] < this->minPoint[0])
        this->minPoint[0] = p[0];
    if (p[1] < this->minPoint[1])
        this->minPoint[1] = p[1];
    if (p[0] > this->maxPoint[0])
        this->maxPoint[0] = p[0];
    if (p[1] > this->maxPoint[1])
        this->maxPoint[1] = p[1];
}

void RectangularBounding::translate(const RSAVector &translation) {
    this->minPoint += translation;
    this->maxPoint += translation;
}

void RectangularBounding::expand(double expansion) {
    this->maxPoint += {{expansion, expansion}};
}

RectangularBounding RectangularBoundingBuilder::forArch(const RSAVector &zeroAngleVector, double angleFrom,
                                                        double angleTo)
{
    Expects(angleFrom >= 0);
    Expects(angleTo >= angleFrom);

    RectangularBounding rectangularBounding;
    rectangularBounding.addPoint(Matrix<2, 2>::rotation(angleFrom) * zeroAngleVector);

    double quarterAngle = M_PI/2;
    while (quarterAngle < angleFrom)
        quarterAngle += M_PI/2;

    recurseQuarters(rectangularBounding, zeroAngleVector, angleTo, quarterAngle);
    return rectangularBounding;
}

RegularDiskopolygonAttributes::RegularDiskopolygonAttributes(const std::string &attr) {
    std::istringstream attrStream(attr);
    std::string type;
    attrStream >> type;
    if (type == "standard") {
    	attrStream >> this->nSides >> this->sideLength >> this->radius;
    } else if (type == "short") {
    	attrStream >> this->nSides >> this->radius;
    	this->sideLength = 2*std::sin(M_PI/this->nSides);
    } else {
        throw ValidationException("Supported argument types are: standard, short");
    }

    ValidateMsg(attrStream, "Malformed attributes. Expected:\n"
                            "standard [number of sites] [side length] [disk radius] 	or\n"
                            "short [number of sides] [disk radius]");
    Validate(this->nSides >= 3);
    Validate(this->sideLength >= 0);
    Validate(this->radius >= 0);

    this->height = this->sideLength / 2 / std::tan(M_PI / this->nSides);
    this->halfDiagonal = this->sideLength / 2 / std::sin(M_PI / this->nSides);
    this->normalizeVolume();
}

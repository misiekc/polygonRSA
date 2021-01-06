//
// Created by Piotr Kubala on 26/03/2020.
//

#ifndef RSA3D_REGULARDISKOPOLYGON_H
#define RSA3D_REGULARDISKOPOLYGON_H

#include "../../Shape.h"
#include "../SpheroCylinder2D.h"

class RectangularBounding {
private:
    RSAVector minPoint{{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()}};
    RSAVector maxPoint{{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}};

public:
    RSAVector getBottomLeft() const { return minPoint; }
    RSAVector getBottomRight() const { return {{maxPoint[0], minPoint[1]}}; }
    RSAVector getTopRight() const { return maxPoint; }
    RSAVector getTopLeft() const { return {{minPoint[0], maxPoint[1]}}; }

    void addPoint(const RSAVector &p);
    void translate(const RSAVector &translation);
    void expand(double expansion);
};

class RectangularBoundingBuilder {
private:
    static void recurseQuarters(RectangularBounding &bounding, const RSAVector &zeroAngleVector, double angleTo,
                                double quarterAngle);

public:
    static RectangularBounding forArch(const RSAVector &zeroAngleVector, double angleFrom, double angleTo);
};

class RegularDiskopolygonAttributes {
private:
    std::size_t nSides{};
    double sideLength{};
    double radius{};
    double height{};
    double halfDiagonal{};

    void normalizeVolume();

public:
    RegularDiskopolygonAttributes() = default;

    /**
     * @brief Calculates attributes from 2 available formats.
     * @details The formats are:
     * <pre>standard [number of sides] [side length] [radius]</pre>
     * <pre>short [number of sides] [radius]</pre>
     * In the second case, the side length is like in the polygon, whose circumscribed circle has radius 1.
     */
    explicit RegularDiskopolygonAttributes(const std::string &attr);

    size_t getNSides() const { return nSides;}
    double getSideLength() const { return sideLength; }
    double getRadius() const { return radius; }
    double getHeight() const { return height; }
    double getHalfDiagonal() const { return halfDiagonal; }
};

class RegularDiskopolygon : public RSAShape {
private:
    static RegularDiskopolygonAttributes attributes;

    SpheroCylinder2D getSpherocylinder(std::size_t index) const;

protected:
    void setOrientation(const RSAOrientation &orientation) override;

public:
    static void initClass(const std::string &attr);

    static std::size_t getNSides() { return attributes.getNSides(); }
    static double getSideLength() { return attributes.getSideLength(); }
    static double getRadius() { return attributes.getRadius(); }
    static double getHeight() { return attributes.getHeight(); }
    static double getHalfDiagonal() { return attributes.getHalfDiagonal(); }

    bool overlap(RSABoundaryConditions *bc, const RSAShape *s) const override;
    bool voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition, const RSAOrientation &orientation,
                     double spatialSize, double angularSize) const override;
    RSAShape *clone() const override;
    std::string toWolfram() const override;
    double getVolume() const override { return 1; }

    double getAngle() const { return this->getOrientation()[0]; };
};


#endif //RSA3D_REGULARDISKOPOLYGON_H

//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_SPHEROCYLINDER2D_H
#define RSA3D_SPHEROCYLINDER2D_H


#include "../AnisotropicShape2D.h"
#include "../../Vector.h"

class SpheroCylinder2D : public AnisotropicShape2D {
private:
    static double radius;
    static double halfDistance;
    static RSAVector centerVector;

    /**
     * @brief A ShapeStaticInfo for SpheroCylinder2D. This one is used, instead of Shape class static info, because
     * spherocylinder in also used in RegularDiskoPolygon, which overrides Shape's static info.
     */
    static ShapeStaticInfo spherocylinderShapeInfo;

    static double pointDistance2(const RSAVector &pos, double angle, const RSAVector &point);

    bool angleInRange(double angle, double rangeStart, double rangeEnd) const;
    double pointDistance2(const RSAVector &p) const;
    bool withinExclusionZone(const RSAVector &pointPos, double angle) const;

protected:
    void setAngle(double angle) override;

public:
    static void initClass(const std::string & attr);

    /**
     * @brief Calculates static values for attr format: [ratio] OR [distance] [radius]
     */
    static void calculateStatic(const std::string &attr);

    static double getRadius() { return radius; }
    static double getDistance() { return 2*halfDistance; }

    double getVolume() const override;

    Shape *clone() const override;

    bool overlap(BoundaryConditions *bc, const Shape *s) const override;
    bool pointInside(BoundaryConditions *bc, const RSAVector &da, double angleFrom, double angleTo) const override;
    void store(std::ostream &f) const override;
    void restore(std::istream &f) override;
    std::string toString() const override;
    std::string toPovray() const override;
    std::string toWolfram() const override;
};


#endif //RSA3D_SPHEROCYLINDER2D_H

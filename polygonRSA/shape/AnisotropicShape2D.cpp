//
// Created by PKua on 07.02.18.
//

#include "AnisotropicShape2D.h"

bool AnisotropicShape2D::pointInside(RSABoundaryConditions *bc, const RSAVector &position,
                                     const RSAOrientation &orientation, double orientationRange) const {
	return this->pointInside(bc, position, orientation[0], orientation[0]+orientationRange);
}

Matrix<2, 2> AnisotropicShape2D::getRotationMatrix() const {
    return Matrix<2, 2>::rotation(this->getAngle());
}

Matrix<2, 2> AnisotropicShape2D::getAntiRotationMatrix() const {
    return Matrix<2, 2>::rotation(-this->getAngle());
}

void AnisotropicShape2D::normalizeAngleRange(double shapeAngle, double *angleFrom, double *angleTo, double interval) {
    while (*angleFrom < shapeAngle) {
        *angleFrom += interval;
        *angleTo += interval;
    }
    while (*angleFrom > shapeAngle + interval) {
        *angleFrom -= interval;
        *angleTo -= interval;
    }
}

double AnisotropicShape2D::normalizeAngle(double angle, double interval) {
    while (angle < 0)
        angle += interval;
    while (angle > interval)
        angle -= interval;
    return angle;
}

void AnisotropicShape2D::setOrientation(const Orientation<1> &orientation) {
    this->setAngle(orientation[0]);
}

void AnisotropicShape2D::setAngle(double angle) {
    RSAConvexShape::setOrientation({angle});
}

double AnisotropicShape2D::getAngle() const{
    return this->getOrientation()[0];
}

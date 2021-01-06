//
// Created by PKua on 07.02.18.
//

#ifndef RSA3D_ANISOTROPICSHAPE2D_H
#define RSA3D_ANISOTROPICSHAPE2D_H


#include "../Matrix.h"
#include "ConvexShape.h"

/**
 * @brief Convenient wrapper around ConvexShape<2, 1>.
 *
 * This class specializes a bit more general ConvexShape<2, 1> class and acts as a base for all anisotropic 2D shapes.
 * Derived classes can now operate on angles, not arrays of angles and automatic normalization is introduced
 * (see setAngle()). Also, it provides some useful subroutines for rotation matrices.
 *
 * Note, that
 * ConvexShape::pointInside(BoundaryConditions*,const Vector<SPATIAL_DIMENSION>&,const Orientation<ANGULAR_DIMENSION>&,double) const
 * and Shape::setOrientation are replaced by more
 * AnisotropicShape2D::pointInside(BoundaryConditions<2>*,const Vector<2>&,double,double) const and
 * AnisotropicShape2D::setAngle. See those methods documentation for more information.
 */
class AnisotropicShape2D : public RSAConvexShape {
private:
    /* Delegated to AnisotropicShape2D::setAngle(double) */
    void setOrientation(const RSAOrientation &orientation) final;
    /* Delegated to AnisotropicShape2D::pointInside(BoundaryConditions*,double*,double,double). */
    bool pointInside(BoundaryConditions *bc, const RSAVector &position, const RSAOrientation &orientation,
                     double orientationRange) const final;

protected:

    /**
     * @brief Returns rotation matrix for - getAngle() angle.
     * @return rotation matrix for - getAngle() angle
     */
    Matrix<2, 2> getAntiRotationMatrix() const;

    /**
     * @brief Returns rotation matrix for getAngle() angle.
     * @return rotation matrix for getAngle() angle
     */
    Matrix<2, 2> getRotationMatrix() const;

    /**
     * @brief Sets new orientation of a shape. It replaces Shape::setOrientation(double*).
     *
     * The angle is normalized using normalizeAngle(double, double) with getVoxelAngularSize() as an interval.
     *
     * Derived shape classes have to override this method when they want to keep track of shape's orientation. One would
     * then typically write:
     * \code
     * void Derived::setAngle(double angle) {
     *     AnisotropicShape2D::setAngle(angle);
     *     // invoke getAngle and for example compute vertices
     * }
     * \endcode
     * @param angle new orientation of a shape
     */
    virtual void setAngle(double angle);

public:
    /**
     * @brief Add/subtracts integer multiple of @a interval to @a angleFrom and @a angleTo so that @a angleFrom is
     * bigger than the @a shapeAngle, however not more than @a interval.
     * @param shapeAngle angle of the shape to normalize with respect to
     * @param angleFrom angle range beginning
     * @param angleTo angle range ending
     * @param interval interval to fit angleFrom into
     */
    static void normalizeAngleRange(double shapeAngle, double *angleFrom, double *angleTo, double interval);

    /**
     * @brief Add/subtracts integer multiple of @a interval to @a angle so that it lies in (0, @a interval) and returns
     * result.
     * @param angle angle to normalize
     * @param interval interval to fit angle into
     * @return normalized angle
     */
    static double normalizeAngle(double angle, double interval);


    /**
     * @brief Returns the orientation of a shape.
     * @return the orientation of a shape
     */
    double getAngle() const;

    /**
     * @brief Checks if a given point @a da is within intersection of all excluded zone for all orientations in
     * (@a angleFrom, @a angleTo) range. It replaces Shape::pointInside(BoundaryConditions*,double*,double*,double).
     * @param bc boundary conditions to take into account
     * @param da position of a point to check
     * @param angleFrom array of beginnings of angle intervals
     * @param angleTo array of lengths of angle intervals
     * @return 0 if point is outside, nonzero number otherwise
     */
    virtual bool pointInside(BoundaryConditions *bc, const RSAVector &da, double angleFrom, double angleTo) const = 0;

};


#endif //RSA3D_ANISOTROPICSHAPE2D_H

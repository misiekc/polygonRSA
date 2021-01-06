/*
 * Shape.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef SHAPE_H_
#define SHAPE_H_


#include <string>
#include <ostream>
#include <istream>
#include <array>

#include "../BoundaryConditions.h"
#include "../Positioned.h"
#include "../RND.h"


class ShapeStaticInfo;

/**
 * @brief A @a 2 dimensional shape providing facilities to generate RSA packings.
 *
 * The class ships with parameters to fine-tune packing generation (@see ShapeStaticInfo), for overlap detection,
 * for storing and restoring. For more information see appropriate methods' documentation.
 *
 * Derived classes should also provide:
 * <ul>
 * <li>a method for initializing parameters of generated shapes from a string, which is non-empty if the class
 * describes the whole family of shapes of a specific kind, for example ellipses (they can have different axes ratio).
 * Signature:
 * @code
 * void (*)(const std::string &attr)
 * @endcode
 * The parameters are then stored in a global fields and are used when generating shapes (@see
 * ShapeStaticInfo::setCreateShapeImpl()). Moreover, <strong>this method is obliged to invoke setShapeStaticInfo()
 * </strong> every time it is called (after shape change or @a attr change), which will override current parameters
 * and info and invalidate previously generated Shapes.</li>
 * </ul>
 *
 * Initialization method will be then hard-coded in ShapeFactory in such a way:
 * @code
 * if (shapeName == "SpecificShape") {
 *     SpecificShape::initClass(shapeAttributes); // use initialization method
 * }
 * @endcode
 *
 * <strong>CALLING INITIALIZATION METHOD IS THE FIRST THING TO DO BEFORE USING A CLASS.</strong>
 *
 * <strong>Creating shapes using ShapeFactory is then the preferred way to create them</strong>. Derived classes do not
 * have to provide public constructors or other ways of creating them.
 *
 * Derived classes are <strong>not required</strong> to handle <strong>shapes of sizes other than from the global
 * state</strong>. All methods have unexpected behavior if initializing method has not been called before, however
 * <strong>one can change parameters or current shape</strong> in any moment by calling initialization method again
 * so <strong>derived classes must be able clean up their previous parameters</strong>. If this happens, all previously
 * generated shapes are considered stale and using them leads to an unexpected behaviour, however all new shapes should
 * work properly.
 *
 * non-zero only if voxelInside is capable of dealing with angle-dependent exclusion zones, ie. when generating
 * saturated RSA packings is supported
 */
class Shape : public RSAPositioned {
public:

    /**
     * @brief A pointer to function for creating shapes - taking RND pointer and returning Shape pointer
     */
    using create_shape_fun_ptr = Shape* (*)(RND *rnd);

protected:

    /**
     * @brief The possible results of early rejection
     */
    enum EarlyRejectionResult {
        /** @brief The collision is present. */
        TRUE,

        /** @brief There is no collision. */
        FALSE,

        /** @brief A compex check is required, nothing has been determined. */
        UNKNOWN
    };

private:

    static ShapeStaticInfo shapeStaticInfo;
    static bool earlyRejectionEnabled;

    RSAOrientation orientation;

    EarlyRejectionResult sphereOnVoxelEarlyRejection(const Vector<2> &relativeSpatialCenter,
                                                     double halfSpatialSize) const;
    EarlyRejectionResult furthestCornerEarlyRejection(const Vector<2> &relativeSpatialCenter,
                                                      double halfSpatialSize) const;

protected:

    /**
     * @brief Sets the static information about the shape after validating it.
     *
     * @param shapeStaticInfo the static information about the shape
     */
    static void setShapeStaticInfo(ShapeStaticInfo shapeStaticInfo);

    /**
     * @brief Translates @a second if it is required by @a bc.
     *
     * Ie. @a second is translated by a vector returned by
     * \code
     * bc->getTranslation(result, this, second);
     * \endcode
     * @param bc boundary conditions to apply
     * @param second shape to be translated according to bc
     */
    virtual void applyBC(BoundaryConditions *bc, Shape *second) const;

    /**
     * @brief Sets shape's new orientation.
     *
     * Derived shape classes have to override this method when they want to keep track of
     * shape's orientation. One would then typically write:
     * \code
     * void Derived::setOrientation(const double *orientation) {
     *     Shape::setOrientation(orientation);
     *     // invoke getOrientation and for example compute vertices
     * }
     * \endcode
     * <strong>rotate(double*) method delegates to this, so no orientation changes will be missed.</strong>
     * @param orientation new shape's orientation
     */
    virtual void setOrientation(const RSAOrientation &orientation);


public:
    /**
     * @brief Number of a shape in a packing.
     */
    int no;

    /**
     * @brief Dimensionless time describing when a shape was added to a packing.
     */
    double time;

    /**
     * @brief Constructs default-oriented shape in the origin of the coordinate system.
     */
    Shape();
    virtual ~Shape() = default;

    /**
     * @brief Returns linear size of a cell in a NeighbourGrid.
     *
     * This size should be as small as possible but big enough to avoid overlapping between shapes having centers in
     * cells that are not neighbours.
     * @return linear size of a cell in a NeighbourGrid
     */
    static double getNeighbourListCellSize();

    /**
     * @brief Returns initial linear size of a (cubic) voxel.
     *
     * This size should be as big as possible but shape with the center inside the voxel have to cover the whole voxel.
     * @return initial linear size of a (cubic) voxel
     */
    static double getVoxelSpatialSize();

    /**
     * @brief Returns angular size of a voxel.
     *
     * Default implementation returns 2 * M_PI. This should be as small as possible but big enough so that angle range
     * (0, getVoxelAngularSize()) describes all possible shape's orientations.
     * @return angular size of a voxel
     */
    static double getAngularVoxelSize();

    /**
     * @brief Returns the radius of the cirbumscribed sphere.
     *
     * @return the radius of the cirbumscribed sphere
     */
    static double getCircumsphereRadius();

    /**
     * @brief Returns the radius of the inscribed sphere.
     *
     * @return the radius of the inscribed sphere
     */
    static double getInsphereRadius();

    /**
     * @brief Returns flags indicating if the shape supports saturated packing generation.
     *
     * Default implementation returns false.
     * @return angular size of a voxel
     */
    static bool getSupportsSaturation();

    /**
     * @brief Returns a pointer to function for creating shapes - taking RND pointer and returning Shape pointer.
     * @return a pointer to function for creating shapes - taking RND pointer and returning Shape pointer
     */
    static create_shape_fun_ptr getCreateShapeImpl();

    /**
     * @brief Return all static information about the shape, including voxel sizes, etc.
     *
     * @return all static information about the shape, including voxel sizes, etc.
     */
    static const ShapeStaticInfo &getShapeStaticInfo() {
        return shapeStaticInfo;
    }

    /**
     * @brief Sets a flag which indicated if the early rejection mechanisms should be on.
     *
     * Early rejection is based on circumsphere and inshpere radiuses and can quickly decice whether shapes do overlap
     * or not in most cases. The flag cat be set to false for testing purposes.
     *
     * @param earlyRejectionEnabled a flag which indicated if the early rejection mechanisms should be on
     */
    static void setEarlyRejectionEnabled(bool earlyRejectionEnabled);

    /**
     * @brief Returns an array of all angles describing shape's orientation.
     * @return array describing shape's orietation
     */
    RSAOrientation getOrientation() const;

    /**
     * @brief Increases all shape's angles by respective values from array @a v.
     * @param v an array of angle deltas
     */
    void rotate(const RSAOrientation &v);

    /**
     * @brief Checks if there is overlap with the shape pointed by @a s.
     *
     * If @a s of `this` has a different size than from a global state is may lead to an unexpected behaviour.
     * @param bc boundary conditions to take into account
     * @param s the second shape
     * @return false if there is no overlap, true otherwise
     */
    virtual bool overlap(BoundaryConditions *bc, const Shape *s) const = 0;


    /**
     * @brief Checks if there is overlap with any shape inside vector @a shapes.
     *
     * If @a s of `this` has a different size than from a global state is may lead to an unexpected behaviour.
     * @param bc boundary conditions to take into account
     * @param shapes vector of shapes
     * @return pointer to overlapping shape or nullptr in no overlap detected
     */
    virtual const Shape* overlap(BoundaryConditions *bc, std::vector<const Shape *> *shapes) const;

    /**
     * @brief Makes a quick check for overlap using inscribed and circumscribed sphere radiuses.
     *
     * Early rejection is based on circumsphere and inshpere radiuses and can quickly decice whether shapes do overlap
     * or not in most cases. The function can be optionally used in ovarlap implementation for a speed gain.
     *
     * @param bc boundary conditions used
     * @param s the second shape to check
     * @param shapeStaticInfo_ shape static info to use here; by default it is a global Shape::shapeStaticInfo, but
     * a different one may be used here
     * @return TRUE or FALSE if result has been determined, UNKNOWN if complex check is required
     */
    EarlyRejectionResult
    overlapEarlyRejection(BoundaryConditions *bc, const Shape *s, ShapeStaticInfo shapeStaticInfo_) const;

    EarlyRejectionResult overlapEarlyRejection(BoundaryConditions *bc, const Shape *s) const;

    /**
     * @brief Returns a volume of the shape. dim is a packing dimension which can be smaller than shape dimension. Then
     * the Volume should be an intersection of shape and packing subvolume.
     *
     * Default implementation returns 1. If derived class supports shapes of volume other than 1 this method should be
     * overriden.
     * @return a volume of the shape
     */
    virtual double getVolume() const;

    /**
     * @brief Checks if whole voxel is inside an exclusion zone of @this shape.
     * Specifically, the procedure returns true any virtual paricle with a center and orientation within the voxel will
     * overlap with @a this shape.
     * Default implementation bases on pointInside() and assumes that the excluded volume of @a this is convex.
     *
     * @param bc boundary conditions to take into account
     * @param voxelPosition position of the "left-bottom" corner of the voxel
     * @param orientation orientation of the "left-bottom" corner of the voxel
     * @param spatialSize spatial size of the voxel
     * @param angularSize angular size of the voxel
     * @return true if the voxel is fully covered by the exclusion zone of @a this shape, false otherwise.
     */
    virtual bool voxelInside(BoundaryConditions *bc, const RSAVector &voxelPosition,
                             const RSAOrientation &orientation, double spatialSize,
                             double angularSize) const = 0;

    /**
     * @brief Makes a quick check for a whole voxel being inside an exclusion zone of @this shape using inscribed and
     * circumscribed sphere radiuses.
     *
     * Early rejection is based on circumsphere and inshpere radiuses and can quickly decice whether shape and voxel do
     * overlap or not in most cases. The function can be optionally used in ovarlap implementation for a speed gain.
     *
     * @param bc boundary conditions to take into account
     * @param voxelPosition position of the "left-bottom" corner of the voxel
     * @param orientation orientation of the "left-bottom" corner of the voxel
     * @param spatialSize spatial size of the voxel
     * @param angularSize angular size of the voxel
     * @return true if the voxel is fully covered by the exclusion zone of @a this shape, false otherwise.
     */
    EarlyRejectionResult voxelInsideEarlyRejection(BoundaryConditions *bc, const RSAVector &voxelPosition,
                                                   const RSAOrientation &orientation, double spatialSize,
                                                   double angularSize) const;

    /**
     * @brief ?? Moves the shape towards given shape s.
     * @param s ??
     * @return ??
     */
    virtual double minDistance(const Shape *s) const;

    /**
     * @brief Returns string representation of the shape.
     *
     * Default implementation returns an empty string.
     * @return string representation of the shape
     */
    virtual std::string toString() const;

    /**
     * @brief Returns povray string representation of the shape.
     *
     * Default implementation returns an empty string.
     * @return povray string representation of the shape
     */
    virtual std::string toPovray() const;

    /**
     * @brief Returns Wolfram Mathematica string representation of the shape.
     *
     * Default implementation returns an empty string.
     * @return Wolfram Mathematica string representation of the shape
     */
    virtual std::string toWolfram() const;

    /**
     * @brief Serializes shape and writes it in binary form to f.
     *
     * Derived classes: Other than for restore(), not only can shapes consistent with a global state be stored, if
     * a class supports them.
     *
     * Information are stored in a binary form in the following order:
     * <ul>
     * <li>@a SPATIAL_DIMENSION parameter as `char`</li>
     * <li>if non-zero, @a ANGULAR_DIMENSION parameter as `char`</li>
     * <li>Shape::no parameter as `int`</li>
     * <li>Shape::time parameter as `double`</li>
     * <li>Positioned::getPosition() as an array of `double`</li>
     * <li>if @a ANGULAR_DIMENSION non-zero, Shape::getOrientation() as an array of `double`</li>
     * </ul>
     *
     * Derived classes can store additional information afterwards.
     * @code
     * void Derived::store(std::ostream &f) const {
     *     Shape::store(f);
     *     // store for example axes ratio
     * }
     * @endcode
     * @param f stream to write serialized shape to
     */
    virtual void store(std::ostream &f) const;

    /**
     * @brief Deserializes shape from f.
     *
     * Derived classes: <strong>Only shapes consistent with a global state can be restored.</strong> std::runtime_error
     * should be thrown otherwise.
     *
     * Format of data is described in restore().
     * @param f stream to write serialized shape to
     * @throws std::runtime_error if dimensions are incompatible
     */
    virtual void restore(std::istream &f);

    /**
     * @brief Returns dynamically allocated exact copy of a shape.
     * @return dynamically allocated exact copy of a shape
     */
    virtual Shape *clone() const = 0;
};


/**
 * @brief A helper class that manages static Shape data, like voxel size, etc.
 *
 * It makes sure that the data is correct and everything has been set. The obligatory parameters to be set:
 * <ul>
 * <li> setCircumsphereRadius()
 * <li> setInsphereRadius()
 * <li> setDefaultCreateShapeImpl()
 * </ul>
 *
 * The parameters:
 * <ul>
 * <li> setNeighbourListCellSize()
 * <li> setVoxelSpatialSize()
 * <li> setExclusionZoneMinSpan()
 * <li> setExclusionZoneMaxSpan()
 * </ul>
 * will be set automatically based on spheres radia. They can be set manually, before or after the previous ones, if
 * needed, but for most cases they should be left as is. For more info check specific setters.
 *
 */
class ShapeStaticInfo {
private:
    using create_shape_fun_ptr = typename Shape::create_shape_fun_ptr;
    static constexpr double NOT_SPECIFIED = std::numeric_limits<double>::infinity();

    double circumsphereRadius = NOT_SPECIFIED;
    double insphereRadius = NOT_SPECIFIED;
    double neighbourListCellSize = NOT_SPECIFIED;
    double voxelSpatialSize = NOT_SPECIFIED;
    double voxelAngularSize = 2*M_PI;
    double exclusionZoneMinSpan = NOT_SPECIFIED;
    double exclusionZoneMaxSpan = NOT_SPECIFIED;
    bool supportsSaturation = false;
    create_shape_fun_ptr createShapeImpl = nullptr;

public:
    /**
     * @brief Returns the radius of the cirbumscribed sphere.
     *
     * @return the radius of the cirbumscribed sphere
     */
    double getCircumsphereRadius() const {
        this->throwIfIncomplete();
        return circumsphereRadius;
    }

    /**
     * @brief Returns the radius of the inscribed sphere.
     *
     * @return the radius of the inscribed sphere
     */
    double getInsphereRadius() const {
        this->throwIfIncomplete();
        return insphereRadius;
    }

    /**
     * @brief Returns the linear size of a neighbour list cell.
     *
     * @return the linear size of a neighbour list cell
     */
    double getNeighbourListCellSize() const {
        this->throwIfIncomplete();
        return neighbourListCellSize;
    }

    /**
     * @brief Returns the initial linear spatial size of a voxel.
     *
     * @return the initial linear spatial size of a voxel
     */
    double getVoxelSpatialSize() const {
        this->throwIfIncomplete();
        return voxelSpatialSize;
    }

    /**
     * @brief Returns the initial angular size of a voxel.
     *
     * @return the initial angular size of a voxel
     */
    double getVoxelAngularSize() const {
        this->throwIfIncomplete();
        return voxelAngularSize;
    }

    /**
     * @brief Returns the smallest possible distance to exclusion zone boundary.
     *
     * @return smallest possible distance to exclusion zone boundary
     */
    double getExclusionZoneMinSpan() const {
        this->throwIfIncomplete();
        return exclusionZoneMinSpan;
    }

    /**
     * @brief Returns the biggest possible distance to exclusion zone boundary.
     *
     * @return bigggest possible distance to exclusion zone boundary
     */
    double getExclusionZoneMaxSpan() const {
        this->throwIfIncomplete();
        return exclusionZoneMaxSpan;
    }

    /**
     * @brief Returns a pointer to function which dynamically allocates the current shape objects.
     *
     * @return a pointer to function which dynamically allocates the current shape objects
     */
    create_shape_fun_ptr getCreateShapeImpl() const {
        this->throwIfIncomplete();
        return createShapeImpl;
    }

    /**
     * @brief Returns a flag which indicates if shape supports saturated packings.
     *
     * @return a flag which indicates if shape supports saturated packings
     */
    bool getSupportsSaturation() const {
        this->throwIfIncomplete();
        return supportsSaturation;
    }

    /**
     * @brief Sets the radius of the circumsbribed sphere.
     *
     * If neighbour list cell size or max exclusion zone span have not been set, they are calculated then automatically
     * based on the radius.
     *
     * @param circumsphereRadius the radius of the circumsbribed sphere
     */
    void setCircumsphereRadius(double circumsphereRadius);

    /**
     * @brief Sets the radius of the inscribed sphere.
     *
     * If voxel spatial size or min exclusion zone span has not been set, it is calculated then automatically based on
     * the radius.
     *
     * @param insphereRadius the radius of the inscribed sphere
     */
    void setInsphereRadius(double insphereRadius);

    /**
     * @brief Sets size of a cell in neighbour grid.
     *
     * The size should be as small as possible, but shapes from not neighbouring cells must not overlap. If it has not
     * been specified before or after setting the circumsphere radius, it was calculated automatically. It should be set
     * manually only in special cases, for example when some rotational freedom is blocked, in essence when not all
     * orientations of a voxel relative to particle are allowed.
     *
     * @param size size of a cell in neighbour grid
     */
    void setNeighbourListCellSize(double neighbourListCellSize);

    /**
     * @brief Sets initial size of a voxel.
     *
     * It indicates the largest possible voxel size that is fully covered by the exclusion zone of the shape placed in
     * it. If it has not been specified before or after setting the insphere radius, it was calculated automatically.
     * It should be set manually only in special cases, for example when some rotational freedom is blocked, in essence
     * when not all orientations of a voxel relative to particle are allowed.
     *
     * @param size size of a voxel
     */
    void setSpatialVoxelSize(double voxelSpatialSize);

    /**
     * @brief Sets initial angular size of a voxel.
     *
     * It indicates the range of angular variable [0, size). By default, the angular size is 2*M_PI.
     * @param size angular size of a voxel
     */
    void setAngularVoxelSize(double voxelAngularSize);

    /**
     * @brief Sets the smallest possible distance to exclusion zone boundary.
     *
     * If it has not been specified before or after setting the insphere radius, it was calculated automatically as
     * the worst case of 2 times the insphere radius. For vast majority of shapes, excluding crazy ones, it should be
     * optimal.
     *
     * @param exclusionZoneMinSpan the smallest possible distance to exclusion zone boundary
     */
    void setExclusionZoneMinSpan(double exclusionZoneMinSpan);

    /**
     * @brief Sets the biggest possible distance to exclusion zone boundary.
     *
     * If it has not been specified before or after setting the insphere radius, it was calculated automatically as
     * the worst case of 2 times the circumsphere radius - the exclusion zone can extend as far as that for example for
     * shapes allowing saturated packings for 2 parallel particles. For shapes not allowing saturated packings, for
     * example cubes, the exclusion zone does not go further than circumsphere radius + insphere radius, so in that case
     * it can be optimized manually.
     *
     * @param exclusionZoneMaxSpan the biggest possible distance to exclusion zone boundary
     */
    void setExclusionZoneMaxSpan(double exclusionZoneMaxSpan);

    /**
     * @brief Sets a flag which inform if the shape supports saturated packing generation.
     *
     * It indicates whether it is possible to get saturated packings of this shape. It is used by PackingGenerator to
     * skip unnecessary voxel analysis.
     * @param true is shape supports saturated packings
     */
    void setSupportsSaturation(bool supportsSaturation);

    /**
     * @brief Sets a function which will be used to create new shapes by ShapeFactory.
     *
     * If @a ANGULAR_DIMENSION is zero, but a Shape is not isotropic, a function should choose random orientation from
     * supplied random number generator, usually with isotropic distribution. All generated shapes must be identical
     * (disregarding orientation).
     *
     * If a function simply returns dynamically allocated default-constructed shape, setDefaultCreateShapeImpl() can be
     * used.
     * @param fptr a pointer to function with signature `Shape* (*)(RND *rnd)`
     */
    void setCreateShapeImpl(create_shape_fun_ptr createShapeImpl);

    /**
     * @brief Sets a function which will be used to create new shapes using default constructor.
     * @tparam ConcreteShape derived Shape class to instantiate
     */
    template<typename ConcreteShape>
    void setDefaultCreateShapeImpl() {
        this->createShapeImpl = []([[maybe_unused]] RND *rnd) {
            return static_cast<Shape *>(new ConcreteShape());
        };
    }

    /**
     * @brief Checks if all fields which don't have default values has been set. If not, throws.
     *
     * This includes circumsphere and insphere radiuses, voxel and neighbour grid spatial sizes and create shape
     * function implementation.
     */
    void throwIfIncomplete() const;
};


#endif /* SHAPE_H_ */
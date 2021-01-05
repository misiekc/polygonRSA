/*
 * Shape.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef SHAPE_H_
#define SHAPE_H_

#include "../BoundaryConditions.h"
#include "../Positioned.h"
#include "../RND.h"
#include "../Utils.h"
#include <string>
#include <ostream>
#include <istream>
#include <array>

/**
 * @brief A @a SPATIAL_DIMENSION dimensional shape providing facilities to generate RSA packings.
 *
 * The class ships with parameters for fine-tuning packing generation (eg initial voxel size), for overlap detection,
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
 * The parameters are then stored in a global state and are used when generating shapes (see setCreateShapeImpl()).
 * Moreover, <strong>this method is obliged to invoke setNeighbourListCellSize(), setVoxelSpatialSize(),
 * setVoxelAngularSize() and setCreateShapeImpl() (or setDefaultCreateShapeImpl())</strong> with appropriate arguments
 * every time it is called (after shape change or @a attr change).</li>
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
 * have to provide public constructors.
 *
 * Derived classes are <strong>not required</strong> to handle <strong>shapes of sizes other than from the global
 * state</strong>. All methods have unexpected behavior if initializing method has not been called before, however
 * <strong>one can change parameters</strong> in any moment by calling initialization method again. If this happens,
 * all previously generated shapes are considered stale and using them leads to an unexpected behaviour, however all new
 * shapes should work properly.
 *
 * @tparam SPATIAL_DIMENSION how many dimensions has the space in which a shape lives
 * @tparam ANGULAR_DIMENSION how many orientational degrees of freedom a shape has. However, this parameter should be
 * non-zero only if voxelInside is capable of dealing with angle-dependent exclusion zones, ie. when generating
 * saturated RSA packings is supported
 */
class RSAShape : public RSAPositioned{
public:
    /**
     * @brief A pointer to function for creating shapes - taking RND pointer and returning Shape pointer
     */
    using create_shape_fun_ptr = RSAShape* (*)(RND *rnd);

private:

    static double voxelSpatialSize;
    static double voxelAngularSize;
    static double neighbourListCellSize;
    static bool supportsSaturation;
    static create_shape_fun_ptr createShapeImpl;

    RSAOrientation orientation;

protected:

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
	virtual void applyBC(RSABoundaryConditions *bc, RSAShape *second) const;

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

    /**
     * @brief Sets initial size of a voxel.
     *
     * Derived Shape classes have to use this method in the initialization method (see class description) to indicate
     * the largest possible voxel size that is fully covered by the exclusion zone of the shape placed in it.
     * @param size size of a voxel
     */
    static void setVoxelSpatialSize(double size);

    /**
     * @brief Sets initial angular size of a voxel.
     *
     * Derived Shape classes have to use this method in the initialization method (see class description) to indicate
     * the range of angular variable [0, size). By default, the angular size is 2*M_PI
     * @param size angular size of a voxel
     */
    static void setVoxelAngularSize(double size);

    /**
     * @brief Sets size of a cell in neighbour grid.
     *
     * Derived Shape classes have to use this method in the initialization method (see class description) to indicate
     * the linear size of the call in them neighbour grid. The size should be as small as possible, but shapes from not
     * neighbouring cells must not overlap.
     * @param size size of a cell in neighbour grid
     */
    static void setNeighbourListCellSize(double size);

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
    static void setCreateShapeImpl(create_shape_fun_ptr fptr);

    /**
     * @brief Sets a function which will be used to create new shapes using default constructor.
     * @tparam SPECIFIC_SHAPE derived Shape class to instantiate
     */
    template<typename SPECIFIC_SHAPE>
    static void setDefaultCreateShapeImpl();


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
    RSAShape();
	virtual ~RSAShape() = default;

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
	static double getVoxelAngularSize();

	/**
     * @brief Returns flags indicating if the shape supports saturated packing generation.
     *
     * Default implementation returns false.
     * @return angular size of a voxel
     */
	static bool getSupportsSaturation();

    /**
     * @brief returns a pointer to function for creating shapes - taking RND pointer and returning Shape pointer.
     * @return a pointer to function for creating shapes - taking RND pointer and returning Shape pointer
     */
    static const create_shape_fun_ptr &getCreateShapeImpl();

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
	virtual bool overlap(RSABoundaryConditions *bc,
                         const RSAShape *s) const = 0;


	/**
	 * @brief Checks if there is overlap with any shape inside vector @a shapes.
	 *
	 * If @a s of `this` has a different size than from a global state is may lead to an unexpected behaviour.
	 * @param bc boundary conditions to take into account
	 * @param shapes vector of shapes
	 * @return pointer to overlapping shape or nullptr in no overlap detected
	 */
	virtual const RSAShape*
				overlap(RSABoundaryConditions *bc,
						std::vector<const RSAShape *> *shapes) const;

	/**
     * @brief Returns a volume of the shape.
     *
     * Default implementation returns 1. If derived class supports shapes of volume other than 1 this method should be overriden.
     * @return a volume of the shape
     */
	virtual double getVolume() const;

    /**
     * @brief Checks if whole voxel is inside an exclusion zone of @this shape.
     * Specifically, the procedure returns true any virtual paricle with a center and orientation within the voxel will overlap with @a this shape
     * Default implementation bases on pointInside() and assumes that the excluded volume of @a this is convex.
     *
     * @param bc boundary conditions to take into account
     * @param voxelPosition position of the "left-bottom" corner of the voxel
     * @param orientation orientation of the "left-bottom" corner of the voxel
     * @param spatialSize spatial size of the voxel
     * @param angularSize angular size of the voxel
     * @return true if the voxel is fully covered by the exclusion zone of @a this shape, false otherwise.
     */
	virtual bool voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition,
                             const RSAOrientation &orientation, double spatialSize,
                             double angularSize) const = 0;

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
    virtual RSAShape *clone() const = 0;
};

template<typename SPECIFIC_SHAPE>
void RSAShape::setDefaultCreateShapeImpl() {
    createShapeImpl = []([[maybe_unused]] RND *rnd) {
        return static_cast<RSAShape *>(new SPECIFIC_SHAPE());
    };
}

#endif /* SHAPE_H_ */

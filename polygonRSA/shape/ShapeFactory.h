/*
 * ShapeFactory.h
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#ifndef SHAPEFACTORY_H_
#define SHAPEFACTORY_H_

#include "Shape.h"
#include <string>


/**
 * @brief A class providing a way to create shapes of a type determined at the runtime.
 */
class ShapeFactory {
private:
    static void initShapeClass0(const std::string &sClass, const std::string &attr);

public:
    /**
     * @brief Creates a shape using a function provided by a specific shape class which was chosen in the last
     * invocation of initShapeClass().
     *
     * Using it before calling initShapeClass() leads to undefined behaviour.
     * @param rnd random number generator to use
     * @return a pointer to a new, dynamically allocated Shape
     */
	static Shape* createShape(RND *rnd);

    /**
     * @brief Initializes a shape class based on passed name and attributes.
     *
     * It invokes specific shape initialization method (see Shape class description) determined by @a sClass and
     * passes @a attr attributes. %Parameters of the shape are then calculated and the shape class provides appropriate
     * function for creating shapes to be used in a factory.
     *
     * A list of available shapes depends on @a RSA_SPATIAL_DIMENSION and @a RSA_ANGULAR_DIMENSION preprocessor
     * constants:
     * <ul>
     * <li> `RSA_SPATIAL_DIMENSION == any && RSA_ANGULAR_DIMENSION == 0`:
     *     <ul>
     *     <li>Sphere</li>
     *     <li>OrientedCuboid</li>
     *     </ul>
     * </li>
     * <li> `RSA_SPATIAL_DIMENSION == 2 && RSA_ANGULAR_DIMENSION == 1`:
     *     <ul>
     *     <li>Ellipse</li>
     *     <li>SpheroCylinder2D</li>
     *     <li>Rectangle</li>
     *     </ul>
     * </li>
     * <li> `RSA_SPATIAL_DIMENSION == 3 && RSA_ANGULAR_DIMENSION == 0`:
     *     <ul>
     *     <li>Cuboid</li>
     *     </ul>
     * </li>
     * </ul>
     * @param sClass name of a shape class
     * @param attr attributes to pass to Shape class initialization method (see appropriate shape class description)
     */
	static void initShapeClass(const std::string &sClass, const std::string &attr);
};


#endif /* SHAPEFACTORY_H_ */

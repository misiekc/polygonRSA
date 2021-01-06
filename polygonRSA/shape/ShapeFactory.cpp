/*
 * ShapeFactory.cpp
 *
 *  Created on: 16.04.2017
 *      Author: ciesla
 */

#include "ShapeFactory.h"
#include "shapes/polygon/Polygon.h"
#include "shapes/polygon/Triangle.h"
#include "shapes/polygon/RegularRoundedPolygonGeometric.h"
#include "shapes/polygon/RegularRoundedPolygonAlgebraic.h"
#include "shapes/polygon/RoundedPolygonAlgebraic.h"

struct UnknownShapeException: public std::domain_error {
	UnknownShapeException() :
			domain_error("") {
	}
};

void ShapeFactory::initShapeClass(const std::string &sClass,
		const std::string &attr) {
	std::ostringstream error;
	error << "[ShapeFactory::initShapeClass] ";
	try {
		ShapeFactory::initShapeClass0(sClass, attr);
		return;
	} catch (UnknownShapeException &e) {
		error << "Unknown shape: " << sClass;
	} catch (ValidationException &e) {
		error << sClass << " attributes:" << std::endl;
		error << "    " << attr << std::endl;
		error << "malformed. Reason:" << std::endl;
		error << "    " << e.what();
	}
	die(error.str());
}

void ShapeFactory::initShapeClass0(const std::string &sClass,
		const std::string &attr) {
	// Shapes for any dimension, without angular dimensions
	if (sClass == "Triangle") {
		Triangle::initClass(attr);
		return;
	} else if (sClass == "Polygon") {
        Polygon::initClass(attr);
        return;
    } else if (sClass == "RegularRoundedPolygonGeometric") {
        RegularRoundedPolygonGeometric::initClass(attr);
        return;
    } else if (sClass == "RegularRoundedPolygonAlgebraic") {
        RegularRoundedPolygonAlgebraic::initClass(attr);
        return;
    } else if (sClass == "RoundedPolygonAlgebraic") {
        RoundedPolygonAlgebraic::initClass(attr);
        return;
    }
	throw UnknownShapeException();
}

Shape *ShapeFactory::createShape(RND *rnd) {
	// Fetch appropriate function from Shape
	auto createShapeImpl = Shape::getCreateShapeImpl();
	return createShapeImpl(rnd);
}

//
// Created by pkua on 10.04.19.
//

#ifndef RSA3D_TRIANGLE_H
#define RSA3D_TRIANGLE_H


#include "Polygon.h"

class Triangle : public Polygon {
private:
    static void parseAttributes(const std::string &args, double &a, double &b, double &c);
    static Vector<2> calculateThirdVertex(double a, double b, double c);

public:
    ~Triangle() override = default;

    static void initClass(const std::string &args);
    static std::string preparePolygonAttributes(const std::string &triangleAttr);
};


#endif //RSA3D_TRIANGLE_H

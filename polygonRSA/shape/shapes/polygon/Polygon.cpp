/*
 * Polygon.cpp
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#include "Polygon.h"
#include "../../../Utils.h"
#include <cmath>
#include <sstream>
#include <optional>


#define VOXEL_INSIDE_FULL_ANGLE_ONLY_VERTICES


std::vector<double> Polygon::vertexR;
std::vector<double> Polygon::vertexTheta;
std::vector<std::pair<size_t, size_t>> Polygon::segments;
std::vector<std::pair<size_t, size_t>> Polygon::helperSegments;

double Polygon::calculateCircumscribedCircleRadius(){
    double result = 0.0;
    for(std::pair<size_t, size_t> segment: Polygon::segments){
        if (result < Polygon::vertexR[segment.first])
            result = Polygon::vertexR[segment.first];
        if (result < Polygon::vertexR[segment.second])
            result = Polygon::vertexR[segment.second];
    }
    return result;
}

double Polygon::calculateInscribedCircleRadius(const RSAVector &origin) {
    double result = std::numeric_limits<double>::infinity();
    for (auto segment : Polygon::segments) {
        RSAVector v3 = Polygon::getStaticVertexPosition(segment.first) - origin;
        RSAVector v4 = Polygon::getStaticVertexPosition(segment.second) - origin;
        result = std::min(result, v3.norm());
        result = std::min(result, v4.norm());

        RSAVector edgeVector = v4 - v3;
        bool projectionInsideEdge = ((v3 * edgeVector) * (v4 * edgeVector) < 0);
        if (projectionInsideEdge) {
            double d = std::abs(v3[0] * v4[1] - v3[1] * v4[0]) / (v4 - v3).norm();
            result = std::min(result, d);
        }
    }
    return result;
}

//calculate the area of the triangle made from the origin, vertex i, and vertex j
double Polygon::getTriangleArea(size_t i, size_t j){
    RSAVector vI = Polygon::getStaticVertexPosition(i);
    RSAVector vJ = Polygon::getStaticVertexPosition(j);

    return 0.5*std::abs(vI[0]*vJ[1] - vI[1]*vJ[0]);
}

/**
 * moves the center of the polygon to (0,0)
 */
void Polygon::centerPolygon() {
    RSAVector origin = Polygon::calculateOptimalOrigin();

    for (size_t i=0; i<Polygon::vertexR.size(); i++){
        RSAVector newV = Polygon::getStaticVertexPosition(i) - origin;
        Polygon::vertexR[i] = newV.norm();
        Polygon::vertexTheta[i] = newV.angle();
    }
}

/**
 * @brief Finds optimal origin for largest insphere.
 * @details Uses pole of inaccessibility algorithm known from geography: https://arxiv.org/pdf/1212.3193.pdf
 */
RSAVector Polygon::calculateOptimalOrigin() {
    double xmin, xmax, ymin, ymax;
    xmin = ymin = std::numeric_limits<double>::infinity();
    xmax = ymax = -std::numeric_limits<double>::infinity();

    for (const auto &vertex : Polygon::getStaticVerticesPositions()) {
        if (vertex[0] < xmin)
            xmin = vertex[0];
        if (vertex[0] > xmax)
            xmax = vertex[0];
        if (vertex[1] < ymin)
            ymin = vertex[1];
        if (vertex[1] > ymax)
            ymax = vertex[1];
    }

    RSAVector origin{};
    double insphereRadius{};
    while (std::pow(xmax - xmin, 2) + std::pow(ymax - ymin, 2) > std::pow(INSPHERE_SEARCH_PRECISION, 2)) {
        bool anyPointInsideSampled{};
        for (std::size_t i{}; i <= INSPHERE_SEARCH_DIVISIONS; i++) {
            RSAVector point;
            point[0] = xmin + (xmax - xmin) * i / INSPHERE_SEARCH_DIVISIONS;
            for (std::size_t j{}; j <= INSPHERE_SEARCH_DIVISIONS; j++) {
                point[1] = ymin + (ymax - ymin) * j / INSPHERE_SEARCH_DIVISIONS;
                if (!Polygon::pointInsidePolygon(point, Polygon::getStaticVerticesPositions()))
                    continue;

                anyPointInsideSampled = true;
                double newRadius = Polygon::calculateInscribedCircleRadius(point);
                if (newRadius > insphereRadius) {
                    insphereRadius = newRadius;
                    origin = point;
                }
            }
        }
        Assert(anyPointInsideSampled);

        double newXSize = (xmax - xmin) / INSPHERE_SEARCH_FACTOR;
        double newYSize = (ymax - ymin) / INSPHERE_SEARCH_FACTOR;
        xmin = origin[0] - newXSize / 2;
        xmax = origin[0] + newXSize / 2;
        ymin = origin[1] - newYSize / 2;
        ymax = origin[1] + newYSize / 2;
    }
    return origin;
}

void Polygon::createStarHelperSegments(){
    Polygon::helperSegments.clear();
    Polygon::vertexR.push_back(0.0);
    Polygon::vertexTheta.push_back(0.0);
    for (size_t i=0; i<Polygon::vertexR.size()-1; i++)
        Polygon::helperSegments.emplace_back(i, Polygon::vertexR.size()-1);
}

/**
 * @param args contains information about polygon
 *
 * Data should be separated by only spaces, and should be:
 * nv [xy/rt] v_01 v_02 ... v_(nv-1),1 v_(nv-1),2 ns si_0 ... si_(ns-1) [starHelperSegments / nhs hsi_01 hsi_02 ... hsi_(nhs-1),1 hsi_(nhs-1),2]
 *
 * xy means cartesian coordinates, and rt means polar coordinates
 * nv - number of vertices
 * v_01 v_02 ... v_(nv-1),1 v_(nv-1),2 - subsequent pair of their coordinates, indexed from 0 to ns-1
 * ns - number of segments - some vertices may be reserved exclusively for helper segments
 * si_0 ... si_(ns-1) - indexes of vertices to create segments from; segments are created cyclicly, meaning their ends
 *      will be (si_0, si_1), (si_1, si_2), ... (si_(ns-1), si_0)
 * starHelperSegments - helper segments will be creating automatically from centre of mass to each vertex
 * nhs - number of helper segments, when being specified manually; can be 0
 * hsi_01 hsi_02 ... hsi_(nhs-1),1 hsi_(nhs-1),2 - subsequent pairs of indexes of vertices to create helper segments
 *      from; the first index from pair is segment beginning, the second is segment end
 *
 * Example format of coordinates
 * 4 xy 1 1 1 -1 -1 -1 -1 1 4 0 1 2 3 starHelperSegments
 * or equivalently
 * 4 rt 1.4142 0.7854 1.4142 2.3562 1.4142 3.9270 1.4142 5.4978 4 0 1 2 3 starHelperSegments
 */
void Polygon::initClass(const std::string &args){
    clearOldData();

    std::istringstream in(args);
    Polygon::parseVertices(in);
    Polygon::parseSegments(in);
    Polygon::centerPolygon();
    Polygon::parseHelperSegments(in);
    Polygon::normalizeVolume(in);

    ShapeStaticInfo shapeInfo;
    shapeInfo.setCircumsphereRadius(Polygon::calculateCircumscribedCircleRadius());
    shapeInfo.setInsphereRadius(Polygon::calculateInscribedCircleRadius());
    shapeInfo.setAngularVoxelSize(2*M_PI);
    shapeInfo.setSupportsSaturation(true);
    shapeInfo.setDefaultCreateShapeImpl<Polygon>();

    RSAShape::setShapeStaticInfo(shapeInfo);
}

void Polygon::normalizeVolume(std::istringstream &in) {
    double area;
    in >> area;
    if (in) {
        Validate(area > 0);
    } else {
        ValidateMsg(Polygon::isPolygonConvex(), "Automatic volume normalization for concave polygons is unsupported");
        area = std::accumulate(segments.begin(), segments.end(), 0.0, [](auto a, auto seg) {
            return a + getTriangleArea(seg.first, seg.second);
        });
    }

    std::for_each(vertexR.begin(), vertexR.end(), [area](auto &vR) { vR /= std::sqrt(area); });
}

bool Polygon::isPolygonConvex() {
    RSAVector v0 = getStaticVertexPosition(0);
    RSAVector v1 = getStaticVertexPosition(1);
    RSAVector prevSegment = v1 - v0;
    std::optional<double> prevSpin = std::nullopt;
    for (std::size_t i = 1; i <= segments.size(); i++) {
        auto segmentIdx = segments[i % segments.size()];
        RSAVector currSegment = getStaticVertexPosition(segmentIdx.second)
                                - getStaticVertexPosition(segmentIdx.first);
        double currSpin = currSegment[0] * prevSegment[1] - currSegment[1] * prevSegment[0];
        if (prevSpin.has_value() && *prevSpin * currSpin < 0)
            return false;
        prevSegment = currSegment;
        prevSpin = currSpin;
    }
    return true;
}

void Polygon::clearOldData() {
    vertexR.clear();
    vertexTheta.clear();
    segments.clear();
    helperSegments.clear();
}

void Polygon::parseVertices(std::istringstream &in) {
    std::size_t n;
    std::string format;

    in >> n;
    in >> format;
    ValidateMsg(n >= 3, "At least 3 vertices should be specified");

    for (size_t i=0; i<n; i++){
        double c1, c2;
        in >> c1;
        in >> c2;

        double r, t;
        if(format == "xy"){
            r = std::sqrt(c1*c1 + c2*c2);
            t = std::atan2(c2, c1);
        }else if (format == "rt") {
            r = c1;
            t = c2;
        }else{
            throw ValidationException("Wrong coordinate format. Use rt or xy");
        }
        Polygon::vertexR.push_back(r);
        Polygon::vertexTheta.push_back(t);
    }
}

void Polygon::parseSegments(std::istringstream &in) {
    std::size_t n;
    in >> n;
    ValidateMsg(n >= 3 && n <= vertexR.size(), "Number of segments should be in [3, numOfVertice] range");

    size_t i1;
    in >> i1;
    ValidateMsg(i1 < Polygon::vertexR.size(), "Wrong vertex number in segment");

    // reading segments
    for(size_t i = 1; i<n; i++){
        size_t i2;
        in >> i2;
        ValidateMsg(i2 != i1, "Two identical adjacent vertices in the boundary");
        ValidateMsg(i2 < Polygon::vertexR.size(), "Wrong vertex number in segment");

        Polygon::segments.emplace_back(i1, i2);
        i1 = i2;
    }
    ValidateMsg(Polygon::segments.back().second != Polygon::segments.front().first,
                "Two identical adjacent vertices in the boundary");
    Polygon::segments.emplace_back(Polygon::segments.back().second, Polygon::segments.front().first);
}

void Polygon::parseHelperSegments(std::istringstream &in) {
    std::size_t n;
    if (in.str().find("starHelperSegments") != std::string::npos) {
        std::string starHelperSegments;
        in >> starHelperSegments;
        ValidateMsg(starHelperSegments == "starHelperSegments", "Malformed Polygon attributes");
        Polygon::createStarHelperSegments();
    } else{
        in >> n;
        // reading helper segments
        for(size_t i = 0; i<n; i++){
            size_t i1, i2;
            in >> i1;
            in >> i2;
            ValidateMsg(i1 < Polygon::vertexR.size() && i2 < Polygon::vertexR.size(), "Wrong vertex number in helper segment");
            Polygon::helperSegments.emplace_back(i1, i2);
        }
    }
}


//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
bool Polygon::lineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
    double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
    double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
    double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
    double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

    return (o1*o2 < 0.0) && (o3*o4 < 0.0);
}

//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
bool Polygon::lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4){
    double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
    double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
    double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
    double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

    if ((o1*o2 < 0.0) && (o3*o4 < 0.0)){
        //the voxel center is intersecting
        double d3 = dx + dtheta*l3;
        double d4 = dx + dtheta*l4;
        double do1 = d3*(std::abs(y2 - y1) + std::abs(x2 - x1));
        double do2 = d4*(std::abs(y2 - y1) + std::abs(x2 - x1));
        double do3 = (d4 + d3)*(std::abs(x1 - x4) + std::abs(y1 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));
        double do4 = (d4 + d3)*(std::abs(x2 - x4) + std::abs(y2 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));

        if (do1 < std::abs(o1) && do2 < std::abs(o2) && do3 < std::abs(o3) && do4 < std::abs(o4))
            return true;
        else
            return false;
    }
    else
        return false;
}

double Polygon::getVolume() const{
    return std::accumulate(Polygon::segments.begin(), Polygon::segments.end(), 0.0, [](auto vol, auto seg) {
        return vol + Polygon::getTriangleArea(seg.first, seg.second);
    });
}


double Polygon::segmentPointDistance2(const RSAVector &s1, const RSAVector &s2, const RSAVector &point){
    double segmentLenght2 = (s2-s1).norm2();
    if (segmentLenght2 == 0.0){
        return (s1-point).norm2();
    }
    // Consider the line extending the segment, parameterized as s1 + t (s2 - s1).
    // We find projection of point p onto the line.
    // It falls where t = [(p-s1) . (s2-s1)] / |s2-s1|^2
    // We clamp t from [0,1] to handle points outside the segment s1 s2.
    double t = std::max(0.0, std::min(1.0, ((point - s1)*(s2 - s1)) / segmentLenght2));
    RSAVector projection = s1 + t * (s2 - s1);  // Projection falls on the segment
    return (point - projection).norm2();
}

bool Polygon::overlapComplexCheck(RSAVector &position, double angle, RSAVector &polposition, double polangle) const{
    // prepare segments set to check for overlapping. Only include segments that are not outside the circumsphere of the other shape.
    std::vector<std::pair<RSAVector, RSAVector>> set, polset;
    std::pair<size_t, size_t> segment;
    std::pair<RSAVector, RSAVector> xySegment;

    double circumsphereRadius2 = std::pow(RSAShape::getCircumsphereRadius(), 2);
    for (size_t i = 0; i < Polygon::segments.size() + Polygon::helperSegments.size(); i++){
        if (i<Polygon::segments.size()){
            segment = Polygon::segments[i];
        }else{
            segment = Polygon::helperSegments[i-Polygon::segments.size()];
        }
        xySegment.first[0]  = polposition[0] + Polygon::vertexR[segment.first] * std::cos(Polygon::vertexTheta[segment.first] + polangle);
        xySegment.first[1]  = polposition[1] + Polygon::vertexR[segment.first] * std::sin(Polygon::vertexTheta[segment.first] + polangle);
        xySegment.second[0] = polposition[0] + Polygon::vertexR[segment.second] * std::cos(Polygon::vertexTheta[segment.second] + polangle);
        xySegment.second[1] = polposition[1] + Polygon::vertexR[segment.second] * std::sin(Polygon::vertexTheta[segment.second] + polangle);

        if (Polygon::segmentPointDistance2(xySegment.first, xySegment.second, position) <= circumsphereRadius2)
            polset.push_back(xySegment);

        xySegment.first[0]  = position[0] + Polygon::vertexR[segment.first] * std::cos(Polygon::vertexTheta[segment.first] + angle);
        xySegment.first[1]  = position[1] + Polygon::vertexR[segment.first] * std::sin(Polygon::vertexTheta[segment.first] + angle);
        xySegment.second[0] = position[0] + Polygon::vertexR[segment.second] * std::cos(Polygon::vertexTheta[segment.second] + angle);
        xySegment.second[1] = position[1] + Polygon::vertexR[segment.second] * std::sin(Polygon::vertexTheta[segment.second] + angle);

        if (Polygon::segmentPointDistance2(xySegment.first, xySegment.second, polposition) <= circumsphereRadius2)
            set.push_back(xySegment);
    }

    for (std::pair<RSAVector, RSAVector> polsegment : polset){

        for (std::pair<RSAVector, RSAVector> segment : set){

            if (Polygon::lineLineIntersect(polsegment.first[0], polsegment.first[1], polsegment.second[0], polsegment.second[1],
                                           segment.first[0], segment.first[1], segment.second[0], segment.second[1]))
                return true;
        }
    }
    return false;
}


bool Polygon::overlap(RSABoundaryConditions *bc, const RSAShape *s) const{
    switch (this->overlapEarlyRejection(bc, s)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    Polygon pol = dynamic_cast<const Polygon&>(*s);
    this->applyBC(bc, &pol);

    RSAVector polposition = pol.getPosition();
    RSAVector position = this->getPosition();

    double angle = this->getOrientation()[0];
    double polangle = pol.getOrientation()[0];

    //complex check
    return this->overlapComplexCheck(position, angle, polposition, polangle);
}

bool Polygon::voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition,
                          const RSAOrientation &voxelOrientation, double spatialSize, double angularSize) const {
    if (voxelOrientation[0] > RSAShape::getAngularVoxelSize())
        return true;

    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, voxelOrientation, spatialSize, angularSize)) {
        case TRUE:      return true;
        case FALSE:     return false;
        case UNKNOWN:   break;
    }

    double halfSpatialSize = 0.5*spatialSize;
    RSAVector voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
    RSAVector spatialCenter = voxelPosition + voxelTranslation + RSAVector{{halfSpatialSize, halfSpatialSize}};

    if (voxelInsideFullAngleCheck(spatialCenter, halfSpatialSize))
        return true;

    double halfAngularSize = 0.5*angularSize;
    double angularCenter = voxelOrientation[0] + halfAngularSize;

    return voxelInsideComplexCheck(spatialCenter, halfSpatialSize, angularCenter, halfAngularSize);
}

bool Polygon::voxelInsideFullAngleCheck(const RSAVector &spatialCenter, double halfSpatialSize) const {
    double pushDistance = RSAShape::getInsphereRadius() - halfSpatialSize * M_SQRT2;
    Assert(pushDistance >= 0);

    // For small polygons full check yields faster generation, while for larger (like 20-gons) checking only vertices
    // is better
#ifdef VOXEL_INSIDE_FULL_ANGLE_ONLY_VERTICES
    return this->pointInsidePushedVertices(spatialCenter, pushDistance);
#else
    return this->pointInsidePushedVertices(spatialCenter, pushDistance)
            || this->pointInsidePushedEdges(spatialCenter, pushDistance)
            || this->pointInsidePolygon(spatialCenter, this->getVerticesPositions());
#endif
}

bool Polygon::pointInsidePushedVertices(const RSAVector &point, double pushDistance) const {
    for (std::size_t i = 0; i < Polygon::vertexR.size(); i++) {
        RSAVector vertex = Polygon::getVertexPosition(i);
        if ((vertex - point).norm2() <= pushDistance*pushDistance)
            return true;
    }

    return false;
}

bool Polygon::pointInsidePushedEdges(const RSAVector &point, double pushDistance) const {
    for (auto segment : segments) {
        RSAVector s1 = getVertexPosition(segment.first);
        RSAVector s2 = getVertexPosition(segment.second);

        RSAVector p1 = point - s1;
        RSAVector p2 = point - s2;
        RSAVector s = s2 - s1;

        if (p1 * s < 0 || p2 * s > 0)
            continue;

        RSAVector N{{-s[1], s[0]}};
        if (std::pow(p1 * N, 2) <= pushDistance * pushDistance * N.norm2())
            return true;
    }
    return false;
}

bool Polygon::pointInsidePolygon(const RSAVector &point, const std::vector<RSAVector> &vertices) {
    std::size_t intersectionsOnLeft = 0;
    for (auto segment : Polygon::segments) {
        RSAVector s1 = vertices[segment.first];
        RSAVector s2 = vertices[segment.second];

        if (s1[1] > s2[1])
            std::swap(s1, s2);

        if (point[1] > s2[1] || point[1] < s1[1])
            continue;

        double xIntersection = (s2[0] * (point[1] - s1[1]) - s1[0] * (point[1] - s2[1])) / (s2[1] - s1[1]);
        if (xIntersection < point[0])
            intersectionsOnLeft++;
    }

    return intersectionsOnLeft % 2 == 1;
}

bool Polygon::voxelInsideComplexCheck(const RSAVector &spatialCenter, double halfSpatialSize, double angularCenter,
                                      double halfAngularSize) const {
    RSAVector position = this->getPosition();
    double angle = getOrientation()[0];
    for (size_t i = 0; i < segments.size() + helperSegments.size(); i++) {
        std::pair<size_t, size_t> segment;
        if (i < segments.size()) {
            segment = segments[i];
        } else {
            segment = helperSegments[i - segments.size()];
        }

        double x1 = position[0] + vertexR[segment.first] * cos(vertexTheta[segment.first] + angle);
        double y1 = position[1] + vertexR[segment.first] * sin(vertexTheta[segment.first] + angle);
        double x2 = position[0] + vertexR[segment.second] * cos(vertexTheta[segment.second] + angle);
        double y2 = position[1] + vertexR[segment.second] * sin(vertexTheta[segment.second] + angle);

        for (size_t j = 0; j < segments.size() + helperSegments.size(); j++) {
            std::pair<size_t, size_t> vsegment;
            if (j < segments.size()) {
                vsegment = segments[j];
            } else {
                vsegment = helperSegments[j - segments.size()];
            }
            double x3 = spatialCenter[0] + vertexR[vsegment.first] * cos(vertexTheta[vsegment.first] + angularCenter);
            double y3 = spatialCenter[1] + vertexR[vsegment.first] * sin(vertexTheta[vsegment.first] + angularCenter);
            double x4 = spatialCenter[0] + vertexR[vsegment.second] * cos(vertexTheta[vsegment.second] + angularCenter);
            double y4 = spatialCenter[1] + vertexR[vsegment.second] * sin(vertexTheta[vsegment.second] + angularCenter);
            if (lineVoxelIntersect(x1, y1, x2, y2, x3, y3, x4, y4, halfSpatialSize, halfAngularSize,
                                   vertexR[vsegment.first], vertexR[vsegment.second])) {
                return true;
            }
        }
    }
    return false;
}

RSAShape *Polygon::clone() const {
    return new Polygon(*this);
}

std::string Polygon::toPovray() const{
    std::stringstream out;
    out.precision(std::numeric_limits< double >::max_digits10);
    out << "  polygon {" << Polygon::segments.size()+1 << ", ";
    this->vertexToPovray(Polygon::segments[0].first, out);
    for (size_t i=0; i < Polygon::segments.size(); i++) {
        this->vertexToPovray(Polygon::segments[i].second, out);
        if (i<Polygon::segments.size()-1)
            out << " ,";
    }
    out << "  texture { finish { ambient 1 diffuse 0 } pigment { color Red} } }" << std::endl;

    return out.str();
}

std::string Polygon::toString() const {
    return this->toWolfram();
}

std::string Polygon::toWolfram() const {
    std::ostringstream out;
    out.precision(std::numeric_limits<double>::max_digits10);
    out << "Polygon[{";
    out << this->getVertexPosition(Polygon::segments[0].first);
    for (std::size_t i = 1; i < segments.size(); i++)
        out << ", " << this->getVertexPosition(Polygon::segments[i].first);
    out << "}]";
    return out.str();
}

RSAVector Polygon::getVertexPosition(std::size_t index) const {
    RSAVector position = this->getPosition();
    double angle = Polygon::vertexTheta[index] + this->getOrientation()[0];
    return RSAVector{{
                             position[0] + Polygon::vertexR[index] * std::cos(angle),
                             position[1] + Polygon::vertexR[index] * std::sin(angle)
                     }};
}

std::vector<RSAVector> Polygon::getVerticesPositions() const {
    std::vector<RSAVector> result(Polygon::segments.size());
    for (std::size_t i{}; i < Polygon::segments.size(); i++)
        result[i] = this->getVertexPosition(Polygon::segments[i].first);
    return result;
}

RSAVector Polygon::getStaticVertexPosition(std::size_t index) {
    return RSAVector{{
                             Polygon::vertexR[index] * std::cos(Polygon::vertexTheta[index]),
                             Polygon::vertexR[index] * std::sin(Polygon::vertexTheta[index])
                     }};
}

std::vector<RSAVector> Polygon::getStaticVerticesPositions() {
    std::vector<RSAVector> result(Polygon::segments.size());
    for (std::size_t i{}; i < Polygon::segments.size(); i++)
        result[i] = Polygon::getStaticVertexPosition(Polygon::segments[i].first);
    return result;
}

void Polygon::vertexToPovray(std::size_t index, std::ostream &out) const {
    RSAVector position = this->getVertexPosition(index);
    out << "< " << position[0] << ", " << position[1] << ", 0.0002>";
}
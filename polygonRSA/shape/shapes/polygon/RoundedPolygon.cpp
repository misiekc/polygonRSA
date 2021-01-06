/*
 * Polygon.cpp
 *
 *  Created on: 16.04.2018
 *      Author: ciesla
 */

#include "RoundedPolygon.h"
#include "../../../Utils.h"
#include <cmath>
#include <sstream>


#define VOXEL_INSIDE_FULL_ANGLE_ONLY_VERTICES


double RoundedPolygon::radius;

double RoundedPolygon::calculateCircumscribedCircleRadius(){
	return Polygon::calculateCircumscribedCircleRadius() + RoundedPolygon::radius;
}

void RoundedPolygon::normalizeVolume(std::istringstream &in){
	double area;
	in >> area;
	if (in) {
        Validate(area > 0);
    } else {
        area = RoundedPolygon::getArea();
    }

	std::for_each(vertexR.begin(), vertexR.end(), [area](auto &vR) { vR /= sqrt(area); });
    RoundedPolygon::radius /= sqrt(area);
}

/**
 * @param args contains information about rounded polygon
 *
 * Data should be separated by only spaces, and should be:
 * radius nv [xy/rt] v_01 v_02 ... v_(nv-1),1 v_(nv-1),2 ns si_0 ... si_(ns-1)
 *
 * xy means cartesian coordinates, and rt means polar coordinates
 * nv - number of vertices
 * v_01 v_02 ... v_(nv-1),1 v_(nv-1),2 - subsequent pair of their coordinates, indexed from 0 to ns-1
 * ns - number of segments - some vertices may be reserved exclusively for helper segments
 * si_0 ... si_(ns-1) - indexes of vertices to create segments from with radiuses; segments are created cyclicly, meaning their ends
 *      will be (si_0, si_1), (si_1, si_2), ... (si_(ns-1), si_0)
 *
 * Example format of coordinates
 * 0.2 4 xy 1 1 1 -1 -1 -1 -1 1 4 0 1 2 3
 * or equivalently
 * 0.2 4 rt 1.4142 0.7854 1.4142 2.3562 1.4142 3.9270 1.4142 5.4978 4 0 1 2 3
 */
void RoundedPolygon::initClass(const std::string &args){
    clearOldData();
	std::istringstream in(args);

	in >> RoundedPolygon::radius;

    Polygon::parseVertices(in);
    Polygon::parseSegments(in);
    Polygon::centerPolygon();

    RoundedPolygon::normalizeVolume(in);

    ShapeStaticInfo shapeInfo;
    shapeInfo.setCircumsphereRadius(Polygon::calculateCircumscribedCircleRadius() + RoundedPolygon::radius);
	shapeInfo.setInsphereRadius(Polygon::calculateInscribedCircleRadius() + RoundedPolygon::radius);
	shapeInfo.setAngularVoxelSize(2*M_PI);
	shapeInfo.setSupportsSaturation(true);
	shapeInfo.setDefaultCreateShapeImpl <RoundedPolygon> ();

	RSAShape::setShapeStaticInfo(shapeInfo);
}

double RoundedPolygon::distance2pq(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q){
    double dx = ((x1 + p*(x2 - x1)) - (x3 + q*(x4 - x3)));
    double dy = ((y1 + p*(y2 - y1)) - (y3 + q*(y4 - y3)));
    return dx*dx + dy*dy;
}

void RoundedPolygon::gradientDistance2pq(std::array<double, 2> &gradient, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double p, double q){
    gradient[0] = 2*(x1 - x2)*((p-1)*x1 - p*x2 + x3*(1-q) + q*x4) + 2*(y1 - y2)*((p-1)*y1 - p*y2 + y3*(1-q) + q*y4);
    gradient[1] = 2*(x3 - x4)*((1-p)*x1 + p*x2 + x3*(q-1) - q*x4) + 2*(y3 - y4)*((1-p)*y1 + p*y2 + y3*(q-1) - q*y4);
}

/*
//returns distance between line segment from point 1 to 2 and line segment from point 3 to 4
double RoundedPolygon::lineLineDistance2(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){

	double p = 0.5, q = 0.5;
    double factor = 0.01;
    double precision = 1.0e-8;
    bool precisionTest = false;
    std::array<double, 2> gradient{};
    do{
        gradientDistance2pq(gradient, x1, y1, x2, y2, x3, y3, x4, y4, p, q);
        p -= factor*gradient[0];
        q -= factor*gradient[1];
        if (p<0) p=0;
        if (p>1) p=1;
        if (q<0) q=0;
        if (q>1) q=1;
        if (
            (std::fabs(gradient[0])<precision && std::fabs(gradient[1])<precision) ||
            (std::fabs(gradient[0])<precision && (q==0 || (q==1)) ) ||
            (std::fabs(gradient[1])<precision && (p==0 || (p==1)) ) ||
            ( (p==0 || p==1) && (q==0 || q==1) )
           )
            precisionTest = true;
    }while(!precisionTest);
    return distance2pq(x1, y1, x2, y2, x3, y3, x4, y4, p, q);
}
*/

//returns distance between line segment from point 1 to 2 and line segment from point 3 to 4 for segments that do not intersect
double RoundedPolygon::lineLineDistance2(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){

	RSAVector p1 = {{x1, y1}};
	RSAVector p2 = {{x2, y2}};
	RSAVector p3 = {{x3, y3}};
	RSAVector p4 = {{x4, y4}};

	double d1 = Polygon::segmentPointDistance2(p3, p4, p1);
	double d2 = Polygon::segmentPointDistance2(p3, p4, p2);
	double d3 = Polygon::segmentPointDistance2(p1, p2, p3);
	double d4 = Polygon::segmentPointDistance2(p1, p2, p4);

	return std::min(std::min(std::min(d1,d2),d3),d4);
}

	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4, but endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
bool RoundedPolygon::lineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4){

	if (Polygon::lineVoxelIntersect(x1, y1, x2, y2, x3, y3, x4, y4, dx, dtheta, l3, l4))
		return true;

	double d = std::sqrt(lineLineDistance2(x1, y1, x2, y2, x3, y3, x4, y4));
	double delta = 2*std::max(l3,l4)*std::sin(dtheta/2) + dx*M_SQRT2;

	if (d+delta < 2*RoundedPolygon::radius)
		return true;
	else
		return false;
}

bool RoundedPolygon::overlapComplexCheck(RSAVector &position, double angle, RSAVector &polposition, double polangle) const{
	// prepare segments set to check for overlapping. Only include segments that are not outside the circumsphere of the other shape.
	std::vector<std::pair<RSAVector, RSAVector>> set, polset;
	std::pair<size_t, size_t> segment;
	std::pair<RSAVector, RSAVector> xySegment;

	double circumsphereRadius2 = std::pow(RSAShape::getCircumsphereRadius() + RoundedPolygon::radius, 2);
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
			else if(RoundedPolygon::lineLineDistance2(polsegment.first[0], polsegment.first[1], polsegment.second[0], polsegment.second[1],
					segment.first[0], segment.first[1], segment.second[0], segment.second[1]) < 4*RoundedPolygon::radius*RoundedPolygon::radius)
				return true;
		}
	}
	return false;
}

bool RoundedPolygon::voxelInsideComplexCheck(const RSAVector &spatialCenter, double halfSpatialSize, double angularCenter, double halfAngularSize) const {
/*
	if (
			(std::fabs(spatialCenter[0] - 0.007926502214) < halfSpatialSize) &&
			(std::fabs(spatialCenter[1] - 7.637121042) < halfSpatialSize) &&
			(std::fabs(angularCenter - 5.152016127) < halfAngularSize)
			)
		std::cout << std::endl;
*/
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

double RoundedPolygon::getArea() {
    ValidateMsg(Polygon::isPolygonConvex(), "Automatic volume normalization for concave polygons is unsupported");

    double area = std::accumulate(segments.begin(), segments.end(), 0.0, [](auto a, auto seg) {
        return a + getTriangleArea(seg.first, seg.second);
    });

    size_t i0, i1, i2;
    for(size_t i=0; i<Polygon::segments.size(); i++){
    	i1 = (i>0) ? Polygon::segments[(i-1)].first : Polygon::segments[Polygon::segments.size()-1].first;
    	i0 = Polygon::segments[i].first;
    	i2 = Polygon::segments[(i+1)%Polygon::segments.size()].first;

    	RSAVector v10 = Polygon::getStaticVertexPosition(i1) - Polygon::getStaticVertexPosition(i0);
        RSAVector v20 = Polygon::getStaticVertexPosition(i2) - Polygon::getStaticVertexPosition(i0);
        RSAVector v12 = Polygon::getStaticVertexPosition(i1) - Polygon::getStaticVertexPosition(i2);

        double areaTmp = 0.5*std::abs(v10[0]*v20[1] - v10[1]*v20[0]);
        double h = 2*areaTmp / v12.norm();
        double angle = std::acos(h/v10.norm()) + std::acos(h/v20.norm());

        area += RoundedPolygon::radius*RoundedPolygon::radius*(M_PI - angle)/2.0;
        area += v10.norm()*RoundedPolygon::radius;
    }
    return area;
}

double RoundedPolygon::getVolume() const {
    return 1;
}

RSAShape *RoundedPolygon::clone() const {
    return new RoundedPolygon(*this);
}

std::string RoundedPolygon::toPovray() const{
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	out << "  polygon {" << Polygon::segments.size()+1 << ", ";
	this->vertexToPovray(Polygon::segments[0].first, out);
	out << ", ";
	for (size_t i=0; i < Polygon::segments.size(); i++) {
        this->vertexToPovray(Polygon::segments[i].second, out);
        if (i<Polygon::segments.size()-1)
        	out << " ,";
    }
	out << "  texture { finish { ambient 1 diffuse 0 } pigment { color Red} } }" << std::endl;

	for (size_t i=0; i < Polygon::segments.size(); i++) {

		RSAVector vSegment, vNormal, vBegin, vEnd;
		vBegin = this->getVertexPosition(Polygon::segments[i].first);
		vEnd = this->getVertexPosition(Polygon::segments[i].second);
		out << "  disc { <" << vBegin[0] << ", " << vBegin[1] << ", 0.0002 >,  <0.0, 0.0, 1.0>, " << RoundedPolygon::radius << " texture { finish { ambient 1 diffuse 0 } pigment { color Red } } }" << std::endl;
		vSegment = vEnd - vBegin;
		vNormal[0] = vSegment[1]; vNormal[1] = -vSegment[0];
		vNormal = vNormal.normalized();
		out << "  polygon { 5, "
				<< "<" << (vBegin[0] + vNormal[0]*RoundedPolygon::radius) << ", " << (vBegin[1] + vNormal[1]*RoundedPolygon::radius) << ", 0.0002 >, "
				<< "<" << (vEnd[0] + vNormal[0]*RoundedPolygon::radius) << ", " << (vEnd[1] + vNormal[1]*RoundedPolygon::radius) << ", 0.0002 >, "
				<< "<" << (vEnd[0] - vNormal[0]*RoundedPolygon::radius) << ", " << (vEnd[1] - vNormal[1]*RoundedPolygon::radius) << ", 0.0002 >, "
				<< "<" << (vBegin[0] - vNormal[0]*RoundedPolygon::radius) << ", " << (vBegin[1] - vNormal[1]*RoundedPolygon::radius) << ", 0.0002 >, "
				<< "<" << (vBegin[0] + vNormal[0]*RoundedPolygon::radius) << ", " << (vBegin[1] + vNormal[1]*RoundedPolygon::radius) << ", 0.0002 > "
				<< " texture { finish { ambient 1 diffuse 0 } pigment { color Red} } }" << std::endl;
    }

	return out.str();
}

std::string RoundedPolygon::toWolfram() const {
    std::ostringstream out;
    out << "{" << Polygon::toWolfram() << "";

    for (std::size_t i{}; i < Polygon::segments.size(); i++) {
        RSAVector beg = this->getVertexPosition(Polygon::segments[i].first);
        RSAVector end = this->getVertexPosition(Polygon::segments[i].second);
        out << ", StadiumShape[" << "{" << beg << ", " << end << "}, " << radius << "]";
    }

    out << "}";
	return out.str();
}


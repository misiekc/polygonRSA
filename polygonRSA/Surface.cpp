/*
 * Surface.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Surface.h"
#include <iostream>

Surface::Surface(double s, double ndx) : RSABoundaryConditions() {
	this->size = s;
	this->list = new NeighbourGrid<const RSAShape>(s, ndx);
}

Surface::~Surface() {
	delete this->list;
}

void Surface::clear(){
	this->list->clear();
}


void Surface::add(const RSAShape *s) {
	this->list->add(s, s->getPosition());
}

const RSAShape* Surface::check(const RSAShape *s){
	std::vector<const RSAShape *> neighbours;
	this->list->getNeighbours(&neighbours, s->getPosition());

	return s->overlap(this, &neighbours);
/*
	for(const RSAShape *shape: neighbours) {
		if (shape->overlap(this, s))
			return shape;
	}
	return nullptr;
*/
}


const RSAShape *Surface::getClosestNeighbour(const RSAVector &da, const std::vector<const RSAShape*> &neighbours){
    double d, dmin = std::numeric_limits<double>::max();
    const RSAShape *pmin = nullptr;
    for(const RSAShape *p : neighbours){
        d = this->distance2(da, p->getPosition());
        if (d<dmin){
            pmin = p;
            dmin = d;
        }
    }
    return pmin;
}


const RSAShape *Surface::getClosestNeighbour(const RSAVector &da) {
    std::vector<const RSAShape*> result;
    this->list->getNeighbours(&result, da);
    return getClosestNeighbour(da, result);
}


void Surface::getNeighbours(std::vector<const RSAShape*> *result, const RSAVector &da) {
	this->list->getNeighbours(result, da);
}

NeighbourGrid<const RSAShape> *Surface::getNeighbourGrid(){
	return this->list;
}

RSAVector Surface::checkPosition(const RSAVector &da) const {
	return da;	// do nothing
}

double Surface::distance2(const RSAVector &a1, const RSAVector &a2) const {
	return this->vector(a1 - a2).norm2();
}

RSAVector Surface::vectorFreeBC(const RSAVector &v) const {
	return v;
}

RSAVector Surface::vectorPeriodicBC(const RSAVector &v) const {
    RSAVector result = v;
	for (int i = 0; i < 2; i++) {
		if (v[i] > this->size / 2.0)
			result[i] -= this->size;
		else if (v[i] < -this->size / 2.0)
			result[i] += this->size;
	}
	return result;
}

/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include "Shape.h"
#include <iostream>

double RSAShape::voxelSpatialSize = 0;
double RSAShape::voxelAngularSize = 2*M_PI;
double RSAShape::neighbourListCellSize = 0;

typename RSAShape::create_shape_fun_ptr
        RSAShape::createShapeImpl = nullptr;


void RSAShape::setVoxelSpatialSize(double size) {
	RSAShape::voxelSpatialSize = size;
}

void RSAShape::setVoxelAngularSize(double size) {
	RSAShape::voxelAngularSize = size;
}

void RSAShape::setNeighbourListCellSize(double size) {
	RSAShape::neighbourListCellSize = size;
}

double RSAShape::getVoxelSpatialSize() {
	return RSAShape::voxelSpatialSize;
}

double RSAShape::getVoxelAngularSize() {
	return RSAShape::voxelAngularSize;
}

double RSAShape::getNeighbourListCellSize() {
	return RSAShape::neighbourListCellSize;
}

const typename RSAShape::create_shape_fun_ptr
RSAShape::getCreateShapeImpl() {
    return createShapeImpl;
}

RSAShape::RSAShape() : RSAPositioned(){
    this->orientation.fill(0);
	this->no = 0;
	this->time = 0.0;
}

const RSAShape *
RSAShape::overlap(RSABoundaryConditions *bc, std::vector<const RSAShape *> *shapes) const{

	for(const RSAShape *s: *shapes)
		if (this->overlap(bc, s))
			return s;

	return nullptr;
}



double RSAShape::getVolume() const {
    return 1.;
}

RSAOrientation RSAShape::getOrientation() const{
    return this->orientation;
}

void RSAShape::setOrientation(
		const RSAOrientation &orientation) {
	this->orientation = orientation;
}

void RSAShape::rotate(const RSAOrientation &v){
    RSAOrientation orientation;
	orientation[0] = this->getOrientation()[0] + v[0];
    this->setOrientation(orientation);
}

double RSAShape::minDistance(const RSAShape *s) const{
	return 0.0;
}

std::string RSAShape::toString() const{
	return std::to_string(this->no);
}

std::string RSAShape::toPovray() const{
	return "";
}

std::string RSAShape::toWolfram() const{
	return "";
}

void RSAShape::store(std::ostream &f) const{
	f.write((char *)(&this->no), sizeof(int));
	f.write((char *)(&this->time), sizeof(double));

	// TODO move to Positioned::store?
    double position[2];
    this->getPosition().copyToArray(position);
	f.write((char *)(position), 2*sizeof(double));
	f.write((char *)(this->orientation.data()), sizeof(double));
}

void RSAShape::restore(std::istream &f){
	f.read((char *)(&this->no), sizeof(int));
	f.read((char *)(&this->time), sizeof(double));

	// TODO move to Positioned::restore?
	double position[2];
	f.read((char *)(position), 2*sizeof(double));
	this->setPosition(RSAVector(position));
	f.read((char *)(this->orientation.data()), sizeof(double));
}

void RSAShape::applyBC(RSABoundaryConditions *bc, RSAShape *second) const {
    second->translate(bc->getTranslation(this->getPosition(), second->getPosition()));
}

void RSAShape::setCreateShapeImpl(RSAShape *(*const fptr)(RND *)) {
    createShapeImpl = fptr;
}

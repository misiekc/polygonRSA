/*
 * Voxel.cpp
 *
 *  Created on: 08.03.2017
 *      Author: ciesla
 */

#include "Voxel.h"
#include "Utils.h"

#include <sstream>
#include <stdexcept>
#include <array>
#include <iostream>


Voxel::Voxel(){
	this->lastAnalyzed = 0;
	this->depth = 0;
}


Voxel::Voxel(const RSAVector &pos, const RSAOrientation &angle){
	this->position = pos;
	this->orientation = angle;
	this->lastAnalyzed = 0;
	this->depth = 0;
}

bool Voxel::isInside(const RSAVector &pos, double size){
	for(int i=0; i<2; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}
	return true;
}


bool Voxel::isInside(const RSAVector &pos, double size, const RSAOrientation &angle,
                     double asize) {
	for(int i=0; i<2; i++){
		if (pos[i]<this->position[i])
			return false;
		if (pos[i]>=(this->position[i]+size))
			return false;
	}

	if (angle[0]<this->orientation[0])
		return false;
	if (angle[0]>=(this->orientation[0]+asize))
		return false;

	return true;
}



const RSAVector &Voxel::getPosition(){
	return this->position;
}


const RSAOrientation &Voxel::getOrientation(){
	return this->orientation;
}


std::string Voxel::toPovray(double ssize){
	std::stringstream out;

	out.precision(std::numeric_limits< double >::max_digits10);

	RSAVector da = this->getPosition();
	double x1 = da[0], x2 = da[0] + ssize;
	double y1 = da[1], y2 = da[1] + ssize;

	out << "polygon {5, < " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>, "
					+ "< " + std::to_string(x1) + ", " + std::to_string(y2) + ", 1.0>, "
					+ "< " + std::to_string(x2) + ", " + std::to_string(y2) + ", 1.0>, "
					+ "< " + std::to_string(x2) + ", " + std::to_string(y1) + ", 1.0>, "
					+ "< " + std::to_string(x1) + ", " + std::to_string(y1) + ", 1.0>"
					+ " texture { pigment { color Black} } }";

	/*
			sRes += "\r\n  text { ttf \"timrom.ttf\" \"" + std::to_string(i) + "\" 1, 0 pigment { color White } scale 0.1 translate < ";
			for(unsigned char j=0; j<this->dimension; j++)
				sRes += std::to_string(da[j]+0.5*this->voxelSize) + ", ";
			sRes +=  "1.0003> }";
	*/

	return out.str();
}



std::string Voxel::toWolfram(double ssize, double asize){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	RSAVector da = this->getPosition();
	double x1 = da[0], x2 = da[0] + ssize;
	double y1 = da[1], y2 = da[1] + ssize;

	out << "Polygon[{ {" << x1 << ", " << y1 << "}, "
			<< "{" << x1 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y2 << "}, "
			<< "{" << x2 << ", " << y1 << "} }]";
	out << "(* angles: [ " << this->getOrientation()[0] << ", " << (this->getOrientation()[0] + asize) << ") *)" << std::endl;
	return out.str();
}



std::string Voxel::toString(){
	std::stringstream out;
	out.precision(std::numeric_limits< double >::max_digits10);

	out << " position: (";
	out << this->position[0] << ", " << this->position[1];
	out << ") orientation: (";
	out << this->orientation[0];
	out << ")";
	return out.str();
}


void Voxel::store(std::ostream &f) const{
	double posArray[2];
	this->position.copyToArray(posArray);
	f.write((char *)(posArray), 2*sizeof(double));
	f.write((char *)(this->orientation.data()), sizeof(double));
}


void Voxel::restore(std::istream &f){
	double posArray[2];
	f.read((char *)(posArray), 2*sizeof(double));
	this->position = RSAVector(posArray);
	f.read((char *)(this->orientation.data()), sizeof(double));
}


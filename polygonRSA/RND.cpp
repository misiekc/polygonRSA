/*
 * RND.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "RND.h"

RND::RND(int seed) {
	this->mt = new std::mt19937(seed);
	this->distribution = new std::uniform_real_distribution<double>(0.0, 1.0);
}

RND::RND() {
	this->mt = new std::mt19937();
	this->distribution = new std::uniform_real_distribution<double>(0.0, 1.0);
}

RND::~RND() {
	delete this->mt;
	delete this->distribution;
}

double RND::nextValue(){
	return (*this->distribution)(*this->mt);
}

double RND::nextValue(std::uniform_real_distribution<double> *distr){

	return (*distr)(*this->mt);
}

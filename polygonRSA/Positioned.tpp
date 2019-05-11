/*
 * Positioned.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include<algorithm>
#include "Utils.h"


template <unsigned short SPATIAL_DIMENSION>
const typename Positioned<SPATIAL_DIMENSION>::Offset Positioned<SPATIAL_DIMENSION>::offset;

template <unsigned short SPATIAL_DIMENSION>
const Vector<SPATIAL_DIMENSION> &Positioned<SPATIAL_DIMENSION>::getPosition() const {
	return this->position;
}

template<unsigned short SPATIAL_DIMENSION>
Positioned<SPATIAL_DIMENSION>::Offset::Offset() {
	std::array<int, SPATIAL_DIMENSION> in{};
	in.fill(0);
	int index = 0;
	do {
		std::copy(in.begin(), in.end(), this->offsets[index].begin());
		index++;
	} while(increment(in.data(), SPATIAL_DIMENSION, 1));
}

template<unsigned short SPATIAL_DIMENSION>
const typename Positioned<SPATIAL_DIMENSION>::Offset::specific_vertex_offset &
Positioned<SPATIAL_DIMENSION>::Offset::operator[](std::size_t vertex) const {
	return this->offsets[vertex];
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::setPosition(const Vector<SPATIAL_DIMENSION> &position) {
	this->position = position;
}

template<unsigned short SPATIAL_DIMENSION>
void Positioned<SPATIAL_DIMENSION>::translate(const Vector<SPATIAL_DIMENSION> &v){
	this->setPosition(this->position + v);
}

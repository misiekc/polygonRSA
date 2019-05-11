/*
 * NeighbourGrid.h
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#ifndef NEIGHBOURGRID_H_
#define NEIGHBOURGRID_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include "BoundaryConditions.h"
#include "Utils.h"


template <typename E>
class NeighbourGrid{

private:
	double linearSize;
	int n;
	double dx;

	// contains vectors of cells (vectors with Positioned* inside)
	std::vector<std::vector<E * > * > lists;
	std::vector<std::vector<int> * > neighbouringCells;

	void init(double size, size_t n){
		this->linearSize = size;
		this->n = n;
		this->dx = size/this->n;
		int length = (int) round(pow(this->n, 2));
		this->lists.reserve(length);
		this->neighbouringCells.reserve(length);

		int in[2];
		double da[2];
		int coords[2];

		for(int i=0; i<length; i++){
			this->lists.push_back(new std::vector<E *>);
			this->neighbouringCells.push_back(new std::vector<int>);
			this->neighbouringCells[i]->reserve((4));

			i2position(da, 2, i, this->dx, this->n);
			in[0] = 0; in[1] = 0;
			coordinates(coords, da, 2, this->linearSize, this->dx, this->n);
			do{
				int iCell = neighbour2i(coords, in, 2, 1, this->n);
				this->neighbouringCells[i]->push_back(iCell);
			}while(increment(in, 2, 2));
			// sort and erase to avoid duplicates - important for small packings
			std::sort( neighbouringCells[i]->begin(), neighbouringCells[i]->end() );
			neighbouringCells[i]->erase( std::unique( neighbouringCells[i]->begin(), neighbouringCells[i]->end() ), neighbouringCells[i]->end() );
		}
	}

public:

	NeighbourGrid(double size, double dx){ // @suppress("Class members should be properly initialized")
		Expects(size > 0);
		Expects(dx > 0);

		size_t n = (size_t)(size/dx);
		if (n == 0)
		    throw std::runtime_error("neighbour grid cell too big");
		this->init(size, n);
	}


	NeighbourGrid(double size, size_t n){ // @suppress("Class members should be properly initialized")
		Expects(size > 0);
		Expects(n > 0);

		this->init(size, n);
	}

	virtual ~NeighbourGrid(){
		for(unsigned int i=0; i<this->lists.size(); i++){
				delete this->lists[i];
				delete this->neighbouringCells[i];
			}
	}

	void add(E* s, const RSAVector &da){
	    double array[2];
	    da.copyToArray(array);
		int i = position2i(array, 2, this->linearSize, this->dx, this->n);
		this->lists[i]->push_back(s);
	}

	void remove(E* s, const RSAVector &da){
	    double array[2];
	    da.copyToArray(array);
		int i = position2i(array, 2, this->linearSize, this->dx, this->n);
		typename std::vector<E *>::iterator it;
		if ( (it = std::find(this->lists[i]->begin(), this->lists[i]->end(), s)) != this->lists[i]->end())
			this->lists[i]->erase(it);
	}

	void clear(){
		for(unsigned int i=0; i<this->lists.size(); i++){
			this->lists[i]->clear();
		}
	}

	std::vector<E*> * getCell(const RSAVector &da){
	    double array[2];
	    da.copyToArray(array);
		int i = position2i(array, 2, this->linearSize, this->dx, this->n);
		return this->lists[i];
	}

	void getNeighbours(std::vector<E *> *result, const RSAVector &da) const{
		result->clear();
		std::vector<E *> *vTmp;

		double array[2];
		da.copyToArray(array);
		int i = position2i(array, 2, this->linearSize, this->dx, this->n);
		for(int iCell : *(this->neighbouringCells[i])){
			vTmp = (this->lists[iCell]);
			result->insert(result->end(), vTmp->begin(), vTmp->end());
		}
	}
};

#endif /* NEIGHBOURGRID_H_ */

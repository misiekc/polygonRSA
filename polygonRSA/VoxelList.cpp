/*
 * VoxelList.cpp
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#include "VoxelList.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/**
 * d - requested initial size of a voxel
 */

double VoxelList::findFloorSize(double d){
	double dRes = 1.0;
	while (dRes > d)
		dRes /= 2;
	while (2*dRes < d)
		dRes *= 2;
	return dRes;
}


double VoxelList::findCeilSize(double d){
	double dRes = 1.0;
	while (dRes < d)
		dRes *= 2;
	while (dRes > d)
		dRes /= 2;
	return 2*dRes;
}



inline size_t VoxelList::findArraySize(double range, double cellSize){
	return (size_t)(range/cellSize) + 1;
}

VoxelList::VoxelList(double packingSpatialSize, double requestedSpatialVoxelSize, double shapeAngularRange, double requestedAngularVoxelSize){
	Expects(packingSpatialSize > 0.0);
	Expects(requestedSpatialVoxelSize > 0.0);
	Expects(shapeAngularRange > 0.0);
	Expects(requestedAngularVoxelSize > 0.0);

    this->spatialRange = packingSpatialSize;
	this->spatialVoxelSize = this->findFloorSize(requestedSpatialVoxelSize);
	this->initialVoxelSize = this->spatialVoxelSize;

	this->angularRange = shapeAngularRange;
	this->angularVoxelSize = this->findCeilSize(requestedAngularVoxelSize);
	this->initialAngularVoxelSize = this->angularVoxelSize;

	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);
	this->disabled = false;

	size_t ns = this->findArraySize(this->spatialRange, this->spatialVoxelSize);
	size_t na = this->findArraySize(this->angularRange, this->angularVoxelSize);
	size_t ss = (size_t)(pow(ns, 2)+0.5);
	size_t sa = (size_t)(pow(na, 1)+0.5);
	this->voxels = new Voxel*[ss*sa];
	this->activeTopLevelVoxels = new bool[ss];
	this->voxelNeighbourGrid = new NeighbourGrid<Voxel>(this->spatialVoxelSize*ns, ns);

	this->initVoxels();
	for(size_t i = 0; i<ss; i++){
		this->activeTopLevelVoxels[i] = true;
	}


	this->spatialVoxelSize *= this->dxFactor;
	this->beginningVoxelNumber = ss*sa;
	this->length = this->beginningVoxelNumber;
	this->rebuildNeighbourGrid();
	}


VoxelList::~VoxelList() {
	delete this->spatialDistribution;

	for(size_t i=0; i<this->length; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;
	delete[] this->activeTopLevelVoxels;
	delete this->voxelNeighbourGrid;
}

void VoxelList::disable(){
	this->disabled = true;
}

void VoxelList::initVoxels(){
	size_t ns = this->findArraySize(this->spatialRange, this->spatialVoxelSize);
	size_t na = this->findArraySize(this->angularRange, this->angularVoxelSize);
	RSAVector position;
	RSAOrientation orientation{};
	std::array<int, 2> ins{};
	std::array<int, 1> ina{};

    ins.fill(0);
	int index = 0;
	do{
		for(unsigned char i=0; i<2; i++){
			position[i] = this->spatialVoxelSize*ins[i]; // position point to the "left bottom" corner of a voxel
		}
		ina.fill(0);
		do{
			orientation[0] = this->angularVoxelSize*ina[0]; // orientation point to the "left bottom" corner of a voxel
			this->voxels[index] = new Voxel(position, orientation);
			index++;
		}while(increment(ina.data(), 1, na-1));
	}while(increment(ins.data(), 2, ns-1));
}


void VoxelList::rebuildNeighbourGrid(){
	this->voxelNeighbourGrid->clear();
	for(size_t i=0; i<this->length; i++){
		this->voxelNeighbourGrid->add(this->voxels[i], this->voxels[i]->getPosition());
	}
}


void VoxelList::getNeighbours(std::vector<Voxel *> *result, const RSAVector &da){
	return this->voxelNeighbourGrid->getNeighbours(result, da);
}

void VoxelList::compactVoxelArray(Voxel **list, int &endIndex){
	int beginIndex = 0;

	while(list[endIndex] == nullptr && endIndex > -1)
		(endIndex)--;

	while (beginIndex<endIndex){
		if(list[beginIndex] == nullptr){
			list[beginIndex] = list[endIndex];
			list[endIndex] = nullptr;
		}
		beginIndex++;
		while(list[endIndex] == nullptr && endIndex > -1)
			endIndex--;
	}
}

int VoxelList::getIndexOfTopLevelVoxel(const RSAVector &da){

	double daArray[2];
	da.copyToArray(daArray);
	int n = (int)(this->spatialRange/this->initialVoxelSize) + 1;
	int index = position2i(daArray, 2, n*this->initialVoxelSize, this->initialVoxelSize, n);
	return index;
}


void VoxelList::removeTopLevelVoxel(Voxel *v){
	if (this->disabled)
		return;
	this->activeTopLevelVoxels[this->getIndexOfTopLevelVoxel(v->getPosition())]=false;
}


void VoxelList::checkTopLevelVoxels(){

	RND rnd;
	RSAVector pos;
	RSAOrientation angle{};

	for(size_t i=0; i<this->length; i++){
		Voxel *v = this->voxels[i];
		this->getRandomPositionAndOrientation(&pos, &angle, v, &rnd);
		int index = this->getIndexOfTopLevelVoxel(pos);
		if (this->getIndexOfTopLevelVoxel(v->getPosition())!=index){
			std::cout << "checkTopVoxels problem" << std::endl;
		}
		for(int j=0; j<10; j++){
			this->getRandomPositionAndOrientation(&pos, &angle, v, &rnd);
			if (this->getIndexOfTopLevelVoxel(pos)!=index){
				std::cout << "checkTopVoxels problem" << std::endl;
			}
		}

	}
}

Voxel *VoxelList::getVoxel(const RSAVector &pos, const RSAOrientation &angle){
	std::vector<Voxel *> *vTmp = this->voxelNeighbourGrid->getCell(pos);
	for(Voxel *v : *vTmp){
		if (v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}

Voxel *VoxelList::findVoxel(Voxel **list, size_t listSize, const RSAVector &pos, const RSAOrientation &angle){
	for(size_t i=0; i<listSize; i++){
		Voxel *v = list[i];
		if (v!=nullptr && v->isInside(pos, this->spatialVoxelSize, angle, this->angularVoxelSize)){
			return v;
		}
	}
	return nullptr;
}


void VoxelList::splitVoxel(Voxel *v, double spatialSize, double angularSize, Voxel **vRes){

	unsigned short spatialLoop = 4;
	unsigned short angularLoop = 2;

    RSAVector position = v->getPosition();
    RSAVector vpos = v->getPosition();
    RSAOrientation orientation = v->getOrientation();
    RSAOrientation vangle = v->getOrientation();

    std::array<int, 2> inpos{};
    std::array<int, 1> inangle{};
    inpos.fill(0);
    inangle.fill(0);

	for(unsigned short i=0; i<spatialLoop; i++){
		position[0] = vpos[0] + inpos[0]*spatialSize;
		position[1] = vpos[1] + inpos[1]*spatialSize;
		inangle.fill(0);
		for(unsigned short j=0; j<angularLoop; j++){
			orientation[0] = vangle[0] + inangle[0]*angularSize;
			vRes[i*angularLoop + j] = new Voxel(position, orientation);
			increment(inangle.data(), 1, (unsigned char)1);
		} // for j
		increment(inpos.data(), 2, (unsigned char)1);
	} // for i
}


bool VoxelList::isVoxelInsidePacking(Voxel *v){
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();
	for(unsigned short i=0; i < 2; i++){
		if (vpos[i] >= this->spatialRange){
			return false;
		}
	}
	if (vangle[0] >= this->angularRange){
		return false;
	}
	return true;
}

/**
 * returns true when the whole voxel is inside an exclusion area of any shape in shapes
 * To determine it the method tires to split voxel up to level of maxDepth
 */

bool VoxelList::isVoxelInsideExclusionZone(Voxel *v, double spatialSize, double angularSize,
                                           std::vector<const Shape *> *shapes, BoundaryConditions *bc,
                                           unsigned short depth){
	// if voxel is outside the packing it is inside exclusion zone
//	if (!this->isVoxelInsidePacking(v))
//		return true;
	// otherwise checking

	bool isInside = false;
	for(const Shape *s : *shapes){
		isInside = s->voxelInside(bc, v->getPosition(), v->getOrientation(), spatialSize, angularSize);
		if (isInside)
			break;
	}

	if (isInside || depth == 0){
		return isInside;
	// if cannot determine that it is inside and depth > 0 split and recursively check children
	}else{
		// dzielimy voxel na mniejsze
		double ss = spatialSize / 2.0;
		double as = angularSize / 2.0;
		int arrayLenght = 8;
		Voxel **aVoxels = new Voxel*[ arrayLenght ];
		this->splitVoxel(v, ss, as, aVoxels);
		bool bRes = true;
		// sprawdzamy kazdy z mniejszych
		for(int i=0; i<arrayLenght; i++){
			// jesli choc jeden z mniejszych jest false zwracamy false
			if (!this->isVoxelInsideExclusionZone(aVoxels[i], ss, as, shapes, bc, depth-1)){
				bRes = false;
				break;
			}
		}
		// w przeciwnym razie zwracamy true;

		for(int i=0; i<arrayLenght; i++){
			delete aVoxels[i];
		}
		delete[] aVoxels;

		return bRes;
	}
}


bool VoxelList::isTopLevelVoxelActive(Voxel *v){

	// checking if initial voxel containing v is active (do not have a shape inside)
	int index = this->getIndexOfTopLevelVoxel(v->getPosition());
	return this->activeTopLevelVoxels[index];
}

bool VoxelList::analyzeVoxel(Voxel *v, NeighbourGrid<const Shape> *nl, BoundaryConditions *bc, double spatialSize, double angularSize, unsigned short depth){
	if (!this->disabled){ // && (depth > v->depth || depth==0) ){

	    if (!isTopLevelVoxelActive(v) || !this->isVoxelInsidePacking(v) )
			return true;

	    std::vector<const Shape*> tmpShapes, shapes;
		    nl->getNeighbours(&tmpShapes, v->getPosition());

		int maxNo = v->lastAnalyzed;
		for(const Shape *s: tmpShapes){
			if (v->lastAnalyzed < s->no || depth > v->depth){
				shapes.push_back(s);
				if(maxNo < s->no)
					maxNo = s->no;
			}
		}

		bool isInside = this->isVoxelInsideExclusionZone(v, spatialSize, angularSize, &shapes, bc, depth);

		v->depth = depth;
		v->lastAnalyzed = maxNo;
		return isInside;
	}
	return false;
}


size_t VoxelList::analyzeVoxels(BoundaryConditions *bc, NeighbourGrid<const Shape> *nl, unsigned short depth) {

	size_t begin = this->length;

	_OMP_PARALLEL_FOR
	for (size_t i = 0; i < this->length; i++) {
		Voxel *v = this->voxels[i];
		bool bRemove = this->analyzeVoxel(v, nl, bc, this->spatialVoxelSize, this->angularVoxelSize, depth);
		if (bRemove){
			this->voxelNeighbourGrid->remove(v, v->getPosition());
			delete v;
			this->voxels[i] = nullptr;
		}
		if (i%10000 == 0){ std::cout << "."; std::cout.flush(); }
	}

	std::cout << " compacting" << std::flush;

	int endIndex = this->length-1;
	this->compactVoxelArray(this->voxels, endIndex);
	this->length = endIndex+1;

	return begin - this->length;
}


bool VoxelList::splitVoxels(size_t maxVoxels, NeighbourGrid<const Shape> *nl, BoundaryConditions *bc){
	if (this->disabled)
		return false;
	size_t voxelsFactor = 8;
	if (voxelsFactor*this->length > maxVoxels){
		return false;
	}

	int newListSize = voxelsFactor*(this->length);
	Voxel** newList = new Voxel*[ newListSize ];
	for(int i=0; i<newListSize; i++){
		newList[i] = nullptr;
	}

	Voxel ***aVoxels = new Voxel**[_OMP_MAXTHREADS];
	for(int i=0; i<_OMP_MAXTHREADS; i++){
		aVoxels[i] = new Voxel*[ voxelsFactor ];
	}

	_OMP_PARALLEL_FOR
	for(size_t i=0; i<this->length; i++){
		if (!this->analyzeVoxel(this->voxels[i], nl, bc, this->spatialVoxelSize, this->angularVoxelSize)){ // dividing only not overlapping voxels
			this->splitVoxel(this->voxels[i], this->spatialVoxelSize/2.0, this->angularVoxelSize/2.0, aVoxels[_OMP_THREAD_ID]);
			for(size_t j=0; j<voxelsFactor; j++){
				Voxel *v = aVoxels[_OMP_THREAD_ID][j];
//				if(this->isVoxelInsidePacking(v) && ( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc) ) ){
				if( nl==nullptr || bc==nullptr || !this->analyzeVoxel(v, nl, bc, this->spatialVoxelSize/2.0, this->angularVoxelSize/2.0) ){
					if(this->voxels[i]->depth > 0){
						v->depth = this->voxels[i]->depth-1;
					}
					newList[i*voxelsFactor + j] = v;
				}else{
					delete aVoxels[_OMP_THREAD_ID][j];
				}
			}
		}
		delete this->voxels[i];
		if (i%10000 == 0){ std::cout << "." << std::flush; }
	}

	for(int i=0; i<_OMP_MAXTHREADS; i++){
		delete[] aVoxels[i];
	}
	delete[] aVoxels;

	delete[] this->voxels;

	int endIndex = newListSize - 1;

	std::cout << " compacting" << std::flush;

	this->compactVoxelArray(newList, endIndex);

	this->spatialVoxelSize = (this->spatialVoxelSize/2.0)*this->dxFactor;
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);

	this->angularVoxelSize = (this->angularVoxelSize/2.0)*this->dxFactor;
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	this->length = endIndex+1;
	this->voxels = newList;
	this->rebuildNeighbourGrid();

//	this->checkTopLevelVoxels();
	return true;
}


Voxel * VoxelList::getRandomVoxel(RND *rnd){
	double d = rnd->nextValue();
	return this->voxels[(int)(d*(this->length))];
}


Voxel * VoxelList::getVoxel(int i){
	return this->voxels[i];
}


void VoxelList::getRandomPositionAndOrientation(RSAVector *position, RSAOrientation *orientation, Voxel *v, RND *rnd){
	RSAVector vpos = v->getPosition();
	RSAOrientation vangle = v->getOrientation();

	(*position)[0] = vpos[0] + rnd->nextValue(this->spatialDistribution);
    (*position)[1] = vpos[1] + rnd->nextValue(this->spatialDistribution);
    (*orientation)[0] = vangle[0] + rnd->nextValue(this->angularDistribution);
}


double VoxelList::getVoxelSize(){
	return this->spatialVoxelSize;
}


double VoxelList::getVoxelAngularSize(){
	return this->angularVoxelSize;
}


Voxel* VoxelList::get(int i){
	return this->voxels[i];
}


size_t VoxelList::getLength() const{
	return this->length;
}


double VoxelList::getVoxelsSurface(){
	double result = 0;
	for(size_t i = 0; i< this->length; i++){
		double s = 1.0;
		Voxel *v = this->voxels[i];
		RSAVector position = v->getPosition();
		for(unsigned short j=0; j<2; j++){
			if (position[j]+this->spatialVoxelSize > this->spatialRange){
				s *= this->spatialRange - position[j];
			}else{
				s *= this->spatialVoxelSize;
			}
		}
	    double a = 1.0;
		RSAOrientation orientation = v->getOrientation();
		if (orientation[0]+this->angularVoxelSize > this->angularRange){
			a *= this->angularRange - orientation[0];
		}else{
			a *= this->angularVoxelSize;
		}
		a /= this->angularRange;
		result += s*a;
	}
	return result;
}


std::string VoxelList::toPovray(){
	std::string sRes = "";

	for(size_t i=0; i<this->length; i++){
		sRes += this->voxels[i]->toPovray(this->spatialVoxelSize) + "\r\n";
	}
	return sRes;
}


std::string VoxelList::toWolfram(){
	std::stringstream out;

	for(size_t i=0; i<this->length; i++){
		out << this->voxels[i]->toWolfram(this->spatialVoxelSize, this->angularVoxelSize);
		if (i!=this->length-1)
			out << ", ";
		out << std::endl;
		out << "(* angles: [ " << this->voxels[i]->getOrientation()[0] << ", " << (this->voxels[i]->getOrientation()[0] + this->angularVoxelSize) << ") *)" << std::endl;
	}
	return out.str();
}


void VoxelList::store(std::ostream &f) const{
	f.write((char *) &this->spatialVoxelSize, sizeof(double));
	f.write((char *) &this->angularVoxelSize, sizeof(double));

	size_t size = this->length;
	f.write((char *)(&size), sizeof(size_t));
	for(size_t i=0; i<size; i++){
		this->voxels[i]->store(f);
	}
}


void VoxelList::restore(std::istream &f){
	f.read((char *) &(this->spatialVoxelSize), sizeof(double));
	delete this->spatialDistribution;
	this->spatialDistribution = new std::uniform_real_distribution<double>(0.0, this->spatialVoxelSize);

	f.read((char *)&this->angularVoxelSize, sizeof(double));
	delete this->angularDistribution;
	this->angularDistribution = new std::uniform_real_distribution<double>(0.0, this->angularVoxelSize);

	this->voxelNeighbourGrid->clear();
	for(size_t i=0; i<this->length; i++){
		delete this->voxels[i];
	}
	delete[] this->voxels;

	size_t size;
	int topIndex;
	f.read((char *)&size, sizeof(size_t));
	this->voxels = new Voxel *[size];
	for(size_t i=0; i<size; i++){
		Voxel *v = new Voxel();
		v->restore(f);
		this->voxels[i] = v;
		topIndex = this->getIndexOfTopLevelVoxel(v->getPosition());
		this->activeTopLevelVoxels[topIndex] = true;
	}
	this->length = size;
	this->rebuildNeighbourGrid();
}

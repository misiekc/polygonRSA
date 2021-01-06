/*
 * PackingGenerator.cpp
 *
 *  Created on: 16.04.2017
 *      Author: Michal Ciesla
 */

#include <memory>
#include <chrono>
#include <cstring>
#include "PackingGenerator.h"
#include "surfaces/NBoxPBC.h"
#include "surfaces/NBoxFBC.h"
#include "shape/ShapeFactory.h"
#include "ThreadLocalRND.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <dirent.h>


#ifdef _OPENMP
#include <omp.h>

#endif


double PackingGenerator::FACTOR_LIMIT = 5.0;


PackingGenerator::PackingGenerator(int seed, const Parameters *params) {
	this->params = *params;
	this->seed = seed;
	RND rnd(this->seed);

	this->spatialSize = this->params.surfaceSize;
	this->angularSize = Shape::getAngularVoxelSize();
	if (this->params.requestedAngularVoxelSize > this->angularSize)
		this->params.requestedAngularVoxelSize = this->angularSize;

	this->voxels = new VoxelList(this->spatialSize, Shape::getVoxelSpatialSize(), this->angularSize, this->params.requestedAngularVoxelSize);

	double gridSize = Shape::getNeighbourListCellSize();

	if (this->params.boundaryConditions == "free")
		this->surface = new NBoxFBC(this->params.surfaceSize, gridSize, Shape::getVoxelSpatialSize());
	else
		this->surface = new NBoxPBC(this->params.surfaceSize, gridSize, Shape::getVoxelSpatialSize());
}


PackingGenerator::~PackingGenerator() {
	delete this->voxels;
	delete this->surface;
}


bool PackingGenerator::isSaturated() {
	return (this->voxels->getLength() == 0);
}


double PackingGenerator::getFactor() {
	return this->surface->getArea() / this->voxels->getVoxelsSurface();
}

bool PackingGenerator::isInside(const RSAVector &position, RSAOrientation &orientation){
	if (position[0]>=this->spatialSize || position[0]<0)
		return false;
	if (position[1]>=this->spatialSize || position[1]<0)
		return false;
	if (orientation[0]>=this->angularSize || orientation[0]<0)
		return false;
	return true;

}


void PackingGenerator::testPacking(const Packing &packing, double maxTime){

	int loop = 1000*_OMP_MAXTHREADS;

	std::cout.precision(std::numeric_limits< double >::max_digits10);

	std::cout << "[" << this->seed << " PackingGenerator::testPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	RND rnd(this->seed);
	Shape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume() / this->surface->getArea();
	delete s;

	this->packing = packing;
	for(const Shape *s : packing)
		this->surface->add(s);
	std::cout << "[" << this->seed << " PackingGenerator::testPacking] " << packing.size() << " shapes restored" << std::endl;



	double t = 0;

    ThreadLocalRND threadRND(rnd);

	while (t<maxTime) {

		std::cout << "\r" << "[" << this->seed << " PackingGenerator::testPacking] t=" << std::setprecision(4) << t/maxTime << " choosing " << loop << " shapes..." << std::flush;

		_OMP_PARALLEL_FOR
		for(int i = 0; i<loop; i++){
			Shape *sVirtual = ShapeFactory::createShape(threadRND.get());
			Voxel *v;
			RSAVector pos;
			RSAOrientation angle{};
			do{
				v = this->voxels->getRandomVoxel(threadRND.get());
				this->voxels->getRandomPositionAndOrientation(&pos, &angle, v, threadRND.get());
			}while(!this->isInside(pos, angle));
			// setting shape position and orientation
			sVirtual->translate(pos);
			sVirtual->rotate(angle);
			// checking if shape overlaps with any shape in the packing
			if (this->surface->check(sVirtual)== nullptr){
				_OMP_CRITICAL(stdout)
				{
					std::cout << std::endl << "\t non overlapping shape found " << std::setprecision(10) << sVirtual->toString() << std::endl << std::flush;
					std::cout << "\t povray: " << std::endl << std::setprecision(10) << sVirtual->toPovray() << std::endl << std::flush;
				}
				RSAVector position = sVirtual->getPosition();
				RSAOrientation orientation = sVirtual->getOrientation();
				double delta = 0.0001;
				orientation[0] -= 0.5*delta;

				const Shape *sCovers = nullptr;
				std::vector<const Shape*> vNeighbours;
				this->surface->getNeighbours(&vNeighbours, position);
				for(const Shape *sTmp : vNeighbours){
					if (sTmp->voxelInside(this->surface, position, orientation, 0.0001, delta)){
						sCovers = sTmp;
						break;
					}
				}
				if (sCovers!= nullptr)
				_OMP_CRITICAL(stdout)
				std::cout << "\t in exclusion zone of " << sCovers->toString() << std::endl;
			}
			delete sVirtual;
		} // parallel for

		t += dt * loop;
	} // while

	std::cout << "[" << seed << " PackingGenerator::testPacking] finished after time " << t << std::endl;
}


void PackingGenerator::createPacking() {

	std::cout.precision(std::numeric_limits< double >::max_digits10);
	std::cout << "[" << this->seed << " PackingGenerator::createPacking] using up to " << _OMP_MAXTHREADS;
	std::cout << " concurrent treads" << std::endl;

	std::size_t checkedAgain = 0;
	std::size_t added = 0;
	std::size_t missCounter = 0;
	unsigned short depthAnalyze = 0;

	RND rnd(this->seed);
	Shape *s = ShapeFactory::createShape(&rnd);
	double dt = s->getVolume() / this->surface->getArea();
	delete s;

	int l = 0;
	double t = 0, factor = 1;
	std::size_t tmpSplit = this->params.split, oldTmpSplit = tmpSplit;
//	int snapshotCounter = 0;

    ThreadLocalRND threadRND(rnd);

	const Shape **sOverlapped = new const Shape*[tmpSplit];
	Shape **sVirtual = new Shape*[tmpSplit];
	Voxel **aVoxels = new Voxel *[tmpSplit];

	while (!this->isSaturated()) {

		std::cout << "[" << this->seed << " PackingGenerator::createPacking] choosing " << tmpSplit << " shapes..." << std::flush;
		factor = this->getFactor();
		factor = (factor < 1.0)?1.0:factor;

		_OMP_PARALLEL_FOR
		for(std::size_t i = 0; i<tmpSplit; i++){
			sVirtual[i] = ShapeFactory::createShape(threadRND.get());
			RSAVector pos;
			RSAOrientation angle{};
			do{
				aVoxels[i] = this->voxels->getRandomVoxel(threadRND.get());
				this->voxels->getRandomPositionAndOrientation(&pos, &angle, aVoxels[i], threadRND.get());
			}while(!this->isInside(pos, angle));
			// setting shape position and orientation
			sVirtual[i]->translate(pos);
			sVirtual[i]->rotate(angle);
			// checking if shape overlaps with any shape in the packing
			sOverlapped[i] = this->surface->check(sVirtual[i]);
		} // parallel for


		std::cout << " processing shapes..." << std::flush;

		checkedAgain = 0;
		added = 0;

		// sequentially processing potentially non overlaping shapes
		for(std::size_t i=0; i<tmpSplit; i++){

			t += factor * dt;

			// if there were no intersecting particles in the packing
			if (sOverlapped[i]==nullptr){
				checkedAgain++;
				// checking overlapping again
				sOverlapped[i] = this->surface->check(sVirtual[i]);

				// if there is still no overlap sVirtula is added to the packing and corresponding voxel is removed
				if(sOverlapped[i]== nullptr){ // if non overlapping

					added++;
					l++;
					sVirtual[i]->no = l;
					sVirtual[i]->time = t;

					// consistency check
					if (aVoxels[i]!=this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation())){
						Voxel *v = this->voxels->getVoxel(aVoxels[i]->getPosition(), aVoxels[i]->getOrientation());
						std::cout << std::endl << "Problem: PackingGenerator - inconsistent voxels positions: [" << aVoxels[i]->toString() << "], [" << v->toString() << "]" << std::endl;
						std::cout << "size: " << this->voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << std::endl;
						std::cout << "shape: " << sVirtual[i]->toString() << std::endl;

					}

					this->surface->add(sVirtual[i]);
					this->packing.addShape(sVirtual[i]);
					this->voxels->removeTopLevelVoxel(aVoxels[i]);
				}  //second overlapping check
			} // first overlapping check
			if (sOverlapped[i]!= nullptr){ // removing overlapped virtual shapes
				delete sVirtual[i];
			}
		}  // for

		std::cout << "done, double checked: " << checkedAgain << " added: " << added << ", time: " << t << ", shapes: " << l << std::endl << std::flush;

		//whether splitting voxels
		if (added == 0) {
			missCounter += tmpSplit;
			size_t v0 = this->voxels->getLength(), v1 = v0;

			std::cout << "[" << this->seed << " PackingGenerator::createPacking] splitting " << v0 << " voxels ";
			std::cout.flush();
//			this->toPovray("snapshot_before_" + std::to_string(snapshotCounter++) + ".pov");

			bool b = voxels->splitVoxels(this->params.maxVoxels, this->surface->getNeighbourGrid(), this->surface);
			if (b){
				v1 = this->voxels->getLength();
//				this->toPovray("snapshot_after_" + std::to_string(snapshotCounter++) + ".pov");
				std::cout << " done. " << this->packing.size() << " shapes, " << v1 << " voxels, new voxel size: " << voxels->getVoxelSize() << ", angular size: " << this->voxels->getVoxelAngularSize() << ", factor: " << this->getFactor() << std::endl;
				missCounter = 0;
			}else {
				std::cout << "skipped, analyzing " << this->voxels->getLength() << " voxels, depth = " << depthAnalyze << " " << std::flush;
				this->voxels->analyzeVoxels(this->surface, this->surface->getNeighbourGrid(), depthAnalyze);
				std::cout << " done: " << this->voxels->getLength() << " voxels remained, factor = " << this->getFactor() << std::endl << std::flush;
				tmpSplit = (int)(1.1 * tmpSplit);
				v1 = this->voxels->getLength();
			}
			// if number of voxels has changed
			if (v1!=v0){
				tmpSplit *= ((double)v1 / (double)v0);
			}else{
				tmpSplit = (int)1.1*tmpSplit + _OMP_MAXTHREADS;
			}

			if (tmpSplit > std::max(this->params.maxVoxels/20, 10*this->params.split))
				tmpSplit = std::max(this->params.maxVoxels/20, 10*this->params.split);
			if(v1<v0 && voxels->getLength()<0.001*this->params.maxVoxels && tmpSplit > 10ul*_OMP_MAXTHREADS)
				tmpSplit /= 10.0;
			if(tmpSplit < 10ul*_OMP_MAXTHREADS)
				tmpSplit = 10ul*_OMP_MAXTHREADS;

			if (!b && (double)(v0-v1)/(double)v0 < 0.1){ // not much voxels removed
				depthAnalyze++;
			}else{
				if (depthAnalyze>0)
					depthAnalyze--;
			}

			if(tmpSplit != oldTmpSplit){

				delete[] sOverlapped;
				delete[] sVirtual;
				delete[] aVoxels;

				sOverlapped = new const Shape*[tmpSplit];
				sVirtual = new Shape*[tmpSplit];
				aVoxels = new Voxel *[tmpSplit];

				oldTmpSplit = tmpSplit;
			}
		}else{
			missCounter = 0;
		}
	} // while

	delete[] sOverlapped;
	delete[] sVirtual;
	delete[] aVoxels;

	std::cout << "[" << seed << " PackingGenerator::createPacking] finished after generating " << l << " shapes" << std::endl;
}

void PackingGenerator::run(){
	this->createPacking();
}


const Packing &PackingGenerator::getPacking(){
	return this->packing;
}


void PackingGenerator::toPovray(const Packing &packing, double size, VoxelList *voxels, const std::string &filename){
	std::ofstream file(filename);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera { orthographic location <" << size / 2 << ", " << size / 2 << ", " << (1.3 * size) << "> look_at  <" << size / 2 << ", " << size / 2 << ",  0> }" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size / 2 << ", " << size / 2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;


	for (const Shape *s : packing) {
		file << s->toPovray();
	}


	file << "}" << std::endl;

	if(voxels!=nullptr){
		file << "#declare voxels=union{" << std::endl;
		file << voxels->toPovray() << std::endl;
		file << "}" << std::endl;
	}
	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;

	if (voxels!=nullptr)
		file << "  object { voxels }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}


void PackingGenerator::toPovray(const std::string &filename){
	PackingGenerator::toPovray(this->packing, this->params.surfaceSize, this->voxels, filename);
}


void PackingGenerator::toWolfram(const RSAVector &da, const std::string &filename){

	std::vector<const Shape*> vShapes;
	this->surface->getNeighbours(&vShapes, da);

	std::vector<Voxel *> vVoxels;
	this->voxels->getNeighbours(&vVoxels, da);

	std::ofstream file(filename);
	file << "Graphics[{Red";

	for (const Shape *s : vShapes) {
		file << ", " << std::endl << s->toWolfram();
	}
	if (!vVoxels.empty()){
		file << ", Black, " << std::endl;
		for(Voxel *v: vVoxels){
			file << v->toWolfram(this->voxels->getVoxelSize(), this->voxels->getVoxelAngularSize());
		}
	}
	file << std::endl << "}]" << std::endl;
	file.close();
}


void PackingGenerator::printRemainingVoxels(const std::string &prefix){
	if (this->voxels->getLength()>20)
		return;
	for(size_t i=0; i<this->voxels->getLength(); i++){
		std::string filename(prefix + "_" + std::to_string(i) + ".nb");
		this->toWolfram(this->voxels->getVoxel(i)->getPosition(), filename);
	}
}


void PackingGenerator::toWolfram(const Packing &packing, VoxelList *voxels, const std::string &filename) {
	std::ofstream file(filename);

    file << "Graphics[{Red";

	for (const Shape *s : packing) {
		file << ", " << std::endl << s->toWolfram();
	}

	if(voxels!=nullptr){
		file << ", Black, " << std::endl << voxels->toWolfram();
	}

	file << std::endl << "}]" << std::endl;

	file.close();
}


void PackingGenerator::toWolfram(const std::string &filename){
    PackingGenerator::toWolfram(this->packing, this->voxels, filename);
}


void PackingGenerator::store(std::ostream &f) const{
	std::size_t size = this->packing.size();
	f.write((char *)(&size), sizeof(int));

	for(const Shape *s: this->packing){
		s->store(f);
	}
	this->voxels->store(f);
}


void PackingGenerator::restore(std::istream &f){
	int size;
	RND rnd;
	f.read((char *)(&size), sizeof(int));
	this->packing.clear();
	this->surface->clear();
	for(int i=0; i<size; i++){
		Shape *s = ShapeFactory::createShape(&rnd);
		s->restore(f);
		this->surface->add(s);
		this->packing.addShape(s);
	}
	this->voxels->restore(f);
}

std::vector<std::string> PackingGenerator::findPackingsInDir(const std::string &dirName) {
    std::string prefix = "packing";
    std::string suffix = ".bin";

    DIR *dir = opendir(dirName.c_str());
    dirent *de;
    std::vector<std::string> filenames;
    while ((de = readdir(dir)) != nullptr) {
        std::string filename = de->d_name;
        if (startsWith(filename, prefix) && endsWith(filename, suffix))
            filenames.push_back(dirName + "/" + filename);
    }
    (void) closedir(dir);
    return filenames;
}

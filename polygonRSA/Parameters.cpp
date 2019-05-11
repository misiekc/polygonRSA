/*
 * Parameters.cpp
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#include "Parameters.h"
#include "Utils.h"
#include "Config.h"
#include <iostream>
#include <fstream>


Parameters::Parameters(std::istream &stream) {
	auto config = Config::parse(stream, '=', true);
	for (const auto &key : config.getKeys()){
		if (key == "maxVoxels")					this->maxVoxels = config.getUnsignedLong(key);
		else if (key == "requestedAngularVoxelSize")	this->requestedAngularVoxelSize = config.getDouble(key);
		else if (key == "split") 						this->split = config.getUnsignedLong(key);
		else if (key == "surfaceVolume") 				this->surfaceSize = pow(config.getDouble(key), 0.5);
		else if (key == "storePackings")	 			this->storePackings = config.getString(key) != "false";
		else if (key == "boundaryConditions") 		    this->boundaryConditions = config.getString(key);
		else if (key == "particleType") 				this->particleType = config.getString(key);
		else if (key == "particleAttributes")			this->particleAttributes = config.getString(key);

		else if (key == "from") 						this->from = config.getUnsignedLong(key);
		else if (key == "collectors") 				    this->collectors = config.getUnsignedLong(key);
		else if (key == "ompThreads") 				    this->ompThreads = config.getInt(key);
		else
		    std::cerr << "[Parameters::Parameters] Warning: unknown parameter " << key << std::endl;
	}

	this->validateData();

	#ifdef _OPENMP
	omp_set_num_threads(this->ompThreads);
	#endif
}

Parameters::Parameters(const std::string &fileName) {
    std::ifstream input(fileName, std::ios::in);
    if (!input)
        throw std::runtime_error("Cannot open configuration file: " + fileName);

    (*this) = Parameters{input};    // Yes, it is possible ;)
}

void Parameters::validateData() {
	Validate(maxVoxels > 0);
	Validate(requestedAngularVoxelSize > 0);
	Validate(from >= 0);
	Validate(collectors > 0);
	Validate(split > 0);
	Validate(!std::isnan(surfaceSize) && surfaceSize > 0);
	Validate(ompThreads > 0);
	Validate(!particleType.empty());
}

Parameters::Parameters() {
    this->validateData();       // Make sure default construction is valid
}

std::string Parameters::getPackingSignature() const {
    std::string sPackingFile;
    char buf[20];
    sprintf(buf, "%.0f", pow(this->surfaceSize, 2));
    std::string size(buf);

    std::string particleAttributes;
    if (this->particleAttributes.length() < 100)
        particleAttributes = replaceAll(this->particleAttributes, " ", "_");
    else
        particleAttributes = replaceAll(this->particleAttributes, " ", "_").substr(0, 100);

    return "packing_" + this->particleType + "_" + particleAttributes + "_" + size;
}

double Parameters::sufraceVolume() const {
    return std::pow(this->surfaceSize, 2);
}

bool Parameters::operator==(const Parameters &rhs) const {
    return
           maxVoxels == rhs.maxVoxels &&
           requestedAngularVoxelSize == rhs.requestedAngularVoxelSize &&
           from == rhs.from &&
           collectors == rhs.collectors &&
           split == rhs.split &&
           surfaceSize == rhs.surfaceSize &&
           storePackings == rhs.storePackings &&
           boundaryConditions == rhs.boundaryConditions &&
           particleType == rhs.particleType &&
           particleAttributes == rhs.particleAttributes &&
           ompThreads == rhs.ompThreads;
}

bool Parameters::operator!=(const Parameters &rhs) const {
    return !(rhs == *this);
}

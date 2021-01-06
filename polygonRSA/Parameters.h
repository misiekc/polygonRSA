/*
 * Parameters.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <limits>
#include <string>
#include <cmath>

#include "Utils.h"

class Parameters {
private:
    void validateData();

public:
	std::size_t maxVoxels = 40000000;
	double requestedAngularVoxelSize = 0.3;
	std::size_t from = 0;
	std::size_t collectors = 1;
	std::size_t split = 10000;
	double surfaceSize = pow(10000.0, 0.5);
	bool storePackings = true;

	std::string boundaryConditions = "periodic";
	std::string particleType = "Triangle";
	std::string particleAttributes = "1.0 1.0 1.0";

	int ompThreads = _OMP_MAXTHREADS;

    bool operator==(const Parameters &rhs) const;

    bool operator!=(const Parameters &rhs) const;

    Parameters();
	explicit Parameters(std::istream &stream);
    explicit Parameters(const std::string &fileName);

    std::string getPackingSignature() const;
    double sufraceVolume() const;
};

#endif /* PARAMETERS_H_ */


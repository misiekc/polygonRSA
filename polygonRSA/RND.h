/*
 * RND.h
 *
 *  Created on: 04.04.2017
 *      Author: ciesla
 */

#ifndef RND_H_
#define RND_H_

#include <random>

class RND {

private:
	std::mt19937 *mt;
	std::uniform_real_distribution<double> *distribution;

public:
	RND();
	RND(int seed);
	virtual ~RND();

	double nextValue();
	double nextValue(std::uniform_real_distribution<double> *distr);
};

#endif /* RND_H_ */

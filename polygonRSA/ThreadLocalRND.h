//
// Created by PKua on 20.07.18.
//

#ifndef RSA3D_THREADLOCALRND_H
#define RSA3D_THREADLOCALRND_H


#include "RND.h"

/**
 * @brief Simple wrapper over RND class "splitting" random generator between multiple OMP threads.
 *
 * It creates such many RND instances as the number of currently used OMP threads and feeds them with random seeds from
 * non-overlapping intervals using initial random number generator provided in a constructor. It is useful in
 * OMP-parallelized for loops where random numbers are needed.
 */
class ThreadLocalRND {
private:
    std::vector<RND> rnds;

public:
    /**
     * @brief Constructor supplying omp-thread-local random number generators with random seeds.
     * @param randomSeeds RND generator to generate seeds
     */
    explicit ThreadLocalRND(RND &randomSeeds);

    /**
     * Returns a pointer to omp-thread-local RND generator.
     * @return a pointer to omp-thread-local RND generator
     */
    RND *get();
};


#endif //RSA3D_THREADLOCALRND_H

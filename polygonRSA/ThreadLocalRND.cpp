//
// Created by PKua on 20.07.18.
//

#include "ThreadLocalRND.h"
#include "Utils.h"

ThreadLocalRND::ThreadLocalRND(RND &randomSeeds) {
    this->rnds.reserve(static_cast<std::size_t>(_OMP_MAXTHREADS));
    for (int i = 0; i < _OMP_MAXTHREADS; i++){
        auto seed = static_cast<int>(1000 * (i + randomSeeds.nextValue()));
        this->rnds.emplace_back(seed);
    }
}

RND *ThreadLocalRND::get() {
    return &(this->rnds[_OMP_THREAD_ID]);
}

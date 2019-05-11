//
// Created by PKua on 28.12.17.
//

#ifndef RSA3D_TIMER_H
#define RSA3D_TIMER_H

#include <chrono>

class Timer
{
    using time_point = std::chrono::system_clock::time_point;

    time_point startTime;
    time_point endTime;

public:
    void start();

    void stop();

    template <typename DUR>
    long count() {
        DUR duration = std::chrono::duration_cast<DUR>(endTime - startTime);
        return duration.count();
    }
};

#endif //RSA3D_TIMER_H

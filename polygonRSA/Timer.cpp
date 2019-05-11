//
// Created by PKua on 28.12.17.
//

#include "Timer.h"
#include <bits/unique_ptr.h>

using namespace std::chrono;

void Timer::start() {
    startTime = std::chrono::_V2::system_clock::now();
}

void Timer::stop() {
    endTime = std::chrono::_V2::system_clock::now();
}

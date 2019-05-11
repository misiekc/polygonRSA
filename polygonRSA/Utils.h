/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#ifndef UTILS_C_
#define UTILS_C_

#include "Vector.h"

#ifdef _OPENMP
    #include <omp.h>
    #define __OMP_STRINGIFY__(x) #x

    #define _OMP_PARALLEL_FOR   _Pragma("omp parallel for")
    #define _OMP_ATOMIC         _Pragma("omp atomic")
    #define _OMP_CRITICAL(x)    _Pragma(__OMP_STRINGIFY__(omp critical(x)))
    #define _OMP_MAXTHREADS     omp_get_max_threads()
    #define _OMP_THREAD_ID      omp_get_thread_num()
#else
    #define _OMP_PARALLEL_FOR
    #define _OMP_ATOMIC
    #define _OMP_CRITICAL(x)
    #define _OMP_MAXTHREADS     1
    #define _OMP_THREAD_ID      0
#endif

// Cpp Core Guidelines-style assertions for design by contract
// https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#i6-prefer-expects-for-expressing-preconditions

// Preconditions check (argument validation)
#define Expects(cond) if (!(cond)) throw std::invalid_argument("Precondition " #cond " failed")

// Postconditions check (results assertion)
#define Ensures(cond) if (!(cond)) throw std::runtime_error("Postcondition " #cond " failed")

// Runtime assertion. Why duplicate assert from cassert? Because we don't want to disable is in release mode and
// be more C++ and throw exception
#define Assert(cond) if (!(cond)) throw std::runtime_error("Assertion " #cond " failed")

// Additional macro for validating things like input from file - wrong input shouldn't be considered as assertion fail,
// because it is not programmer's fault ;)

struct ValidationException : public std::domain_error {
    explicit ValidationException(const std::string &msg) : std::domain_error{msg} { }
};

#define Validate(cond) if (!(cond)) throw ValidationException("Validation " #cond " failed")
#define ValidateMsg(cond, msg) if (!(cond)) throw ValidationException(msg)

template <std::size_t N>
using Orientation = std::array<double, N>;

using RSAVector = Vector<2>;
using RSAOrientation = Orientation<1>;

bool increment(int* in, int inlength, int max);
int position2i(const double* da, int dalength, double size, double dx, int n);
void i2position(double* da, int dalength, int index, double dx, int n);
void coordinates(int* result, const double* da, int dalength, double size, double dx, int n);
int neighbour2i(int* coordinates, int* neighbour, int clength, int offset, int n);

// trim from start
std::string &ltrim(std::string &s);
// trim from end
std::string &rtrim(std::string &s);
// trim from both ends
std::string &trim(std::string &s);
// replaces all occurences of search in source by replace
std::string replaceAll(std::string source, const std::string& search, const std::string& replace);
bool endsWith(const std::string& str, const std::string& suffix);
bool startsWith(const std::string& str, const std::string& prefix);
int lastIndexOf(const std::string &s, char target);

void die(const std::string & reason);

double getAngleToOrigin(const Vector<2> &point);
void rotate2D(double* point, double alpha);

inline double P1(double x) {
    return x;
}

inline double P2(double x) {
    return 0.5*(3*x*x - 1);
}

inline double P3(double x) {
    return 0.5*(5*x*x*x - 3*x);
}

inline double P4(double x) {
    return 0.125*(x*x*(35*x*x-30)+3);
}

#endif /* UTILS_C_ */

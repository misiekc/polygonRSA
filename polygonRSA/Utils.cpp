/*
 * Utils.c
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */


#include <cstring>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include "Vector.h"

bool increment(int* in, int inlength, int max){
	if(inlength==0)
		return false;
	int i=0;
	while(i<inlength){
		in[i]++;
		if (in[i]>max && i<inlength-1){
			in[i] = 0;
			i++;
		}else{
			break;
		}
	}
	if(in[inlength-1]>max)
		return false;
	return true;
}


int position2i(const double* da, int dalength, double size, double dx, int n){
	int result = 0;
	int ix;
	for(int i=dalength-1; i>=0; i--){
		if (da[i]<0)
			ix = (int)((da[i] + size)/dx);
		else if (da[i]>=size)
			ix = (int)((da[i] - size)/dx);
		else
			ix = (int)(da[i]/dx);

		if ( (ix+1)*dx == da[i])
			ix++;
		result = n*result + ix;
	}
	return result;
}

void i2position(double* da, int dalength, int index, double dx, int n){
	for(int i=0; i<dalength; i++){
		int iTmp = index % n;
		da[i] = (iTmp+0.5)*dx;
		index /= n;
	}
}

void coordinates(int* result, const double* da, int dalength, double size, double dx, int n){
	for(int i=dalength-1; i>=0; i--){
		if (da[i]<0)
			result[i] = (int)((da[i] + size)/dx);
		else if (da[i]>=size)
			result[i] = (int)((da[i] - size)/dx);
		else
			result[i] = (int)(da[i]/dx);
	}
}

int neighbour2i(int* coordinates, int* neighbour, int clength, int offset, int n){
	int result = 0;

	for(int i=clength-1; i>=0; i--){
		int ix = coordinates[i] + neighbour[i] - offset;
		if (ix>=n)
			ix -= n;
		else if (ix<0)
			ix += n;
		result = n*result + ix;
	}
	return result;
}

/*
	public static int neighbour2i(int index, int[] neighbour, int offset, int n){
		int[] coordinates = new int[neighbour.length];
		int result = 0;
		for(int i=0; i<coordinates.length; i++){
			coordinates[i] = index % n;
			index /= n;

			coordinates[i] += neighbour[i] - offset;
			if (coordinates[i]>=n)
				coordinates[i] -= n;
			else if (coordinates[i]<0)
				coordinates[i] += n;

		}

		for(int i=coordinates.length-1; i>=0; i--){
			result = n*result + coordinates[i];
		}
		return result;
	}
*/

/*
void testda(double* da, int dalength){
	int in[2]; int inlength = 2;
	double dx = 1.5971099293779696;
	int n = 198;
	double size = 316.22776601683796;
	int radius = 1;
	int index = position2i(da, dalength, size, dx, n);
	int coords[2]; int clength = 2;
	coordinates(coords, da, dalength, size, dx, n );

	for(int i=0; i<inlength; i++){
		in[i] = 0;
	}
	do{
//		int i = Commons.neighbour2i(index, in, 1, n);
		int i = neighbour2i(coords, in, clength, 1, n);
		std::cout << i << " ";
	}while(increment(in, inlength, 2*radius));
}


void test(){
	double da1[2] = {35.972999837894264, 298.6595567936803};
	double da2[2] = {35.972999837894264, 298.3801383300324};
	testda(da1, 2);
	std::cout << std::endl;
	testda(da2, 2);
}
*/

// trim from start
std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}

bool endsWith(const std::string &str, const std::string &suffix) {
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

bool startsWith(const std::string &str, const std::string &prefix) {
    return str.size() >= prefix.size() && 0 == str.compare(0, prefix.size(), prefix);
}


std::string replaceAll(std::string source, const std::string& search, const std::string& replace) {
    size_t pos = 0;
    while ((pos = source.find(search, pos)) != std::string::npos) {
         source.replace(pos, search.length(), replace);
         pos += replace.length();
    }
    return source;
}

int lastIndexOf(const std::string &s, char target){
	int ret = -1;
	int curIdx = 0;
	while(s[curIdx] != '\0'){
		if (s[curIdx] == target) ret = curIdx;
		curIdx++;
	}
	return ret;
}

double getAngleToOrigin(const Vector<2> & point) {
	double angle = atan2(point[1], point[0]);
	if(angle < 0)
		return angle + 2 * M_PI;
	else
		return angle;
}

void rotate2D(double* point, double alpha) {
	double cosa = cos(alpha);
	double sina = sin(alpha);
	double x = point[0];
	double y = point[1];
	point[0] = x*cosa - y*sina;
	point[1] = x*sina + y*cosa;
}

void die(const std::string & reason)
{
    std::cerr << reason << std::endl;
    exit(EXIT_FAILURE);
}
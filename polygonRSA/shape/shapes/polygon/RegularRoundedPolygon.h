//
// Created by pkua on 06.04.2020.
//

#ifndef RSA3D_REGULARROUNDEDPOLYGON_H
#define RSA3D_REGULARROUNDEDPOLYGON_H

#include "RoundedPolygon.h"

class RegularRoundedPolygon : public RoundedPolygon {
public:
    /**
     * @param args contains information about rounded polygon. See RegularDiskopolygonAttributes for format.
     */
    static void initClass(const std::string &attr);
};


#endif //RSA3D_REGULARROUNDEDPOLYGON_H

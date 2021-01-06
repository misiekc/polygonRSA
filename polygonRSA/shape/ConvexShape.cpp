//
// Created by PKua on 12.06.18.
//

#include "ConvexShape.h"

bool RSAConvexShape::pointInside(RSABoundaryConditions *bc, const RSAVector &position) const {
    RSAOrientation zeroAngle;
    zeroAngle.fill(0);
    return this->pointInside(bc, position, zeroAngle, 2*M_PI);
}

bool RSAConvexShape::voxelInside(RSABoundaryConditions *bc, const RSAVector &voxelPosition,
                                 const RSAOrientation &orientation, double spatialSize, double angularSize) const
{
    switch(this->voxelInsideEarlyRejection(bc, voxelPosition, orientation, spatialSize, angularSize)) {
        case EarlyRejectionResult::TRUE:      return true;
        case EarlyRejectionResult::FALSE:     return false;
        case EarlyRejectionResult::UNKNOWN:   break;
    }

    int counterSize = 1 << 2;
    bool isInside = true;
    for(int i=0; i<counterSize; i++){
        RSAVector position;
        for(unsigned short j=0; j<2; j++){
            position[j] = voxelPosition[j] + Positioned<2>::offset[i][j]*spatialSize;
        }
        if(!this->pointInside(bc, position, orientation, angularSize)){
            isInside = false;
            break;
        }
    }

    return isInside;
}

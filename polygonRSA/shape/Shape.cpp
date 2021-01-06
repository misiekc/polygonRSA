/*
 * Shape.cpp
 *
 *  Created on: 07.03.2017
 *      Author: Michal Ciesla
 */

#include <iostream>

#include "Shape.h"


ShapeStaticInfo RSAShape::shapeStaticInfo;

bool RSAShape::earlyRejectionEnabled = true;


RSAShape::RSAShape() : RSAPositioned() {
    this->orientation.fill(0);
    this->no = 0;
    this->time = 0.0;
}

const RSAShape *
RSAShape::overlap(RSABoundaryConditions *bc, std::vector<const RSAShape *> *shapes) const {

    for(const RSAShape *s: *shapes)
        if (this->overlap(bc, s))
            return s;

    return nullptr;
}

double RSAShape::getVolume() const {
    return 1.;
}

RSAOrientation RSAShape::getOrientation() const{
    return this->orientation;
}

void RSAShape::setOrientation(const RSAOrientation &orientation) {
    this->orientation = orientation;
}

void RSAShape::rotate(const RSAOrientation &v){
    RSAOrientation orientation;
    for(unsigned short i=0; i<1; i++)
        orientation[i] = this->getOrientation()[i] + v[i];
    this->setOrientation(orientation);
}

double RSAShape::minDistance([[maybe_unused]] const RSAShape *s) const{
    return 0.0;
}

std::string RSAShape::toString() const{
    return std::to_string(this->no);
}

std::string RSAShape::toPovray() const{
    return "";
}

std::string RSAShape::toWolfram() const{
    return "";
}

void RSAShape::store(std::ostream &f) const{
    unsigned short sd = 2;
    unsigned short ad = 1;
    f.write((char *)(&sd), sizeof(unsigned char));
    if (ad>0)
        f.write((char *)(&ad), sizeof(unsigned char));
    f.write((char *)(&this->no), sizeof(int));
    f.write((char *)(&this->time), sizeof(double));
    // TODO move to Positioned::store?
    double position[2];
    this->getPosition().copyToArray(position);
    f.write((char *)(position), 2*sizeof(double));
    if (ad>0)
        f.write((char *)(this->orientation.data()), 1*sizeof(double));
}

void RSAShape::restore(std::istream &f){
    unsigned char sd = 2;
    unsigned char ad = 1;

    f.read((char *)(&sd), sizeof(unsigned char));
    if (f.gcount()==0){ // end of file
        return;
    }
    if (ad > 0)
        f.read((char *)(&ad), sizeof(unsigned char));

    if (sd!=2 || ad!=1){
        throw std::runtime_error(
                std::string("[ERROR] cannot restore Shape: incompatible dimensions: read ")
                + std::to_string(f.gcount())
                + " bytes.");
    }
    f.read((char *)(&this->no), sizeof(int));
    f.read((char *)(&this->time), sizeof(double));

    // TODO move to Positioned::restore?
    double position[2];
    f.read((char *)(position), sd*sizeof(double));
    this->setPosition(Vector<2>(position));

    if (ad>0)
        f.read((char *)(this->orientation.data()), ad*sizeof(double));
}

void RSAShape::applyBC(RSABoundaryConditions *bc, RSAShape *second) const {
    second->translate(bc->getTranslation(this->getPosition(), second->getPosition()));
}

void RSAShape::setShapeStaticInfo(ShapeStaticInfo shapeStaticInfo) {
    shapeStaticInfo.throwIfIncomplete();
    RSAShape::shapeStaticInfo = std::move(shapeStaticInfo);
}

RSAShape::EarlyRejectionResult
RSAShape::overlapEarlyRejection(RSABoundaryConditions *bc, const RSAShape *s, ShapeStaticInfo shapeStaticInfo_) const
{
    if (!earlyRejectionEnabled)
        return UNKNOWN;

    double distance2 = bc->distance2(this->getPosition(), s->getPosition());

    if (distance2 < 4*std::pow(shapeStaticInfo_.getInsphereRadius(), 2))
        return TRUE;
    else if (distance2 > 4*std::pow(shapeStaticInfo_.getCircumsphereRadius(), 2))
        return FALSE;
    else
        return UNKNOWN;
}

RSAShape::EarlyRejectionResult RSAShape::overlapEarlyRejection(RSABoundaryConditions *bc, const RSAShape *s) const {
    return this->overlapEarlyRejection(bc, s, shapeStaticInfo);
}

void RSAShape::setEarlyRejectionEnabled(bool earlyRejectionEnabled) {
    RSAShape::earlyRejectionEnabled = earlyRejectionEnabled;
}

RSAShape::EarlyRejectionResult
RSAShape::voxelInsideEarlyRejection(RSABoundaryConditions *bc, const Vector<2> &voxelPosition,
                                    [[maybe_unused]] const RSAOrientation &orientation, double spatialSize,
                                    [[maybe_unused]] double angularSize) const
{
    if (!earlyRejectionEnabled)
        return UNKNOWN;

    double halfSpatialSize = 0.5*spatialSize;
    auto voxelTranslation = bc->getTranslation(this->getPosition(), voxelPosition);
    auto spatialCenter = voxelPosition + voxelTranslation + Vector<2>(halfSpatialSize);
    auto relativeSpatialCenter = spatialCenter - this->getPosition();

    /* Both early rejection mechanisms */
    /* EarlyRejectionResult result = sphereOnVoxelEarlyRejection(relativeSpatialCenter, halfSpatialSize);
    if (result != UNKNOWN)
        return result;
    else
        return furthestCornerEarlyRejection(relativeSpatialCenter, halfSpatialSize);*/

    /* Sphere on voxel early rejection only */
    /* return sphereOnVoxelEarlyRejection(relativeSpatialCenter, halfSpatialSize); */

    /* Furthest corner early rejection only - seems the best option */
    return furthestCornerEarlyRejection(relativeSpatialCenter, halfSpatialSize);
}

RSAShape::EarlyRejectionResult
RSAShape::sphereOnVoxelEarlyRejection(const Vector<2> &relativeSpatialCenter, double halfSpatialSize) const {
    double spatialCenterDistance2 = relativeSpatialCenter.norm2();
    double voxelCircumsphereRadius = halfSpatialSize * std::sqrt(2);
    double voxelInsphereRadius = halfSpatialSize;

    if (spatialCenterDistance2 < std::pow(shapeStaticInfo.getExclusionZoneMinSpan() - voxelCircumsphereRadius, 2))
        return TRUE;
    else if (spatialCenterDistance2 > std::pow(shapeStaticInfo.getExclusionZoneMaxSpan() - voxelInsphereRadius, 2))
        return FALSE;
    else
        return UNKNOWN;
}

RSAShape::EarlyRejectionResult
RSAShape::furthestCornerEarlyRejection(const Vector<2> &relativeSpatialCenter, double halfSpatialSize) const {
    auto furthestVoxelCorner = relativeSpatialCenter;
    for (std::size_t j = 0; j < 2; j++) {
        if (relativeSpatialCenter[j] > 0)
            furthestVoxelCorner[j] += halfSpatialSize;
        else
            furthestVoxelCorner[j] -= halfSpatialSize;
    }

    double furthestVoxelCornerDistance2 = furthestVoxelCorner.norm2();
    if (furthestVoxelCornerDistance2 < std::pow(shapeStaticInfo.getExclusionZoneMinSpan(), 2))
        return TRUE;
    else if (furthestVoxelCornerDistance2 > std::pow(shapeStaticInfo.getExclusionZoneMaxSpan(), 2))
        return FALSE;
    else
        return UNKNOWN;
}

double RSAShape::getNeighbourListCellSize() {
    return shapeStaticInfo.getNeighbourListCellSize();
}

double RSAShape::getVoxelSpatialSize() {
    return shapeStaticInfo.getVoxelSpatialSize();
}

double RSAShape::getAngularVoxelSize() {
    return shapeStaticInfo.getVoxelAngularSize();
}

double RSAShape::getCircumsphereRadius() {
    return shapeStaticInfo.getCircumsphereRadius();
}

double RSAShape::getInsphereRadius() {
    return shapeStaticInfo.getInsphereRadius();
}

bool RSAShape::getSupportsSaturation() {
    return shapeStaticInfo.getSupportsSaturation();
}

RSAShape::create_shape_fun_ptr RSAShape::getCreateShapeImpl() {
    return shapeStaticInfo.getCreateShapeImpl();
}

void ShapeStaticInfo::setCircumsphereRadius(double circumsphereRadius) {
    Expects(circumsphereRadius > 0);
    Expects(circumsphereRadius < std::numeric_limits<double>::infinity());
    this->circumsphereRadius = circumsphereRadius;
    if (this->neighbourListCellSize == NOT_SPECIFIED)
        this->neighbourListCellSize = 2*circumsphereRadius;
    if (this->exclusionZoneMaxSpan == NOT_SPECIFIED)
        this->exclusionZoneMaxSpan = 2*circumsphereRadius;
}

void ShapeStaticInfo::setInsphereRadius(double insphereRadius) {
    Expects(insphereRadius > 0);
    Expects(insphereRadius < std::numeric_limits<double>::infinity());
    this->insphereRadius = insphereRadius;
    if (this->voxelSpatialSize == NOT_SPECIFIED)
        this->voxelSpatialSize = 2*insphereRadius/std::sqrt(2);
    if (this->exclusionZoneMinSpan == NOT_SPECIFIED)
        this->exclusionZoneMinSpan = 2*insphereRadius;
}

void ShapeStaticInfo::setNeighbourListCellSize(double neighbourListCellSize) {
    Expects(neighbourListCellSize > 0);
    Expects(neighbourListCellSize < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::neighbourListCellSize = neighbourListCellSize;
}

void ShapeStaticInfo::setSpatialVoxelSize(double voxelSpatialSize) {
    Expects(voxelSpatialSize > 0);
    Expects(voxelSpatialSize < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::voxelSpatialSize = voxelSpatialSize;
}

void ShapeStaticInfo::setAngularVoxelSize(double voxelAngularSize) {
    Expects(voxelAngularSize > 0);
    Expects(voxelAngularSize <= 2*M_PI);
    ShapeStaticInfo::voxelAngularSize = voxelAngularSize;
}

void ShapeStaticInfo::setExclusionZoneMinSpan(double exclusionZoneMinSpan) {
    Expects(exclusionZoneMinSpan > 0);
    Expects(exclusionZoneMinSpan < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::exclusionZoneMinSpan = exclusionZoneMinSpan;
}

void ShapeStaticInfo::setExclusionZoneMaxSpan(double exclusionZoneMaxSpan) {
    Expects(exclusionZoneMaxSpan > 0);
    Expects(exclusionZoneMaxSpan < std::numeric_limits<double>::infinity());
    ShapeStaticInfo::exclusionZoneMaxSpan = exclusionZoneMaxSpan;
}

void ShapeStaticInfo::setSupportsSaturation(bool supportsSaturation) {
    ShapeStaticInfo::supportsSaturation = supportsSaturation;
}

void ShapeStaticInfo::setCreateShapeImpl(
        ShapeStaticInfo::create_shape_fun_ptr createShapeImpl) {
    ShapeStaticInfo::createShapeImpl = createShapeImpl;
}

void ShapeStaticInfo::throwIfIncomplete() const {
    if (this->circumsphereRadius == NOT_SPECIFIED)
        throw std::runtime_error("Circumsphere radius has not been set!");
    else if (this->insphereRadius == NOT_SPECIFIED)
        throw std::runtime_error("Insphere radius has not been set!");
    else if (this->createShapeImpl == nullptr)
        throw std::runtime_error("Create shape function implementation radius has not been set!");
}
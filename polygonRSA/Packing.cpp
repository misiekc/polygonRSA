//
// Created by PKua on 17.06.18.
//

#include "Packing.h"
#include "shape/ShapeFactory.h"
#include "Utils.h"
#include <fstream>

Packing::~Packing() {
    for (auto shape : this->packing)
        delete shape;
}

Packing::Packing(const Packing &other) {
    for (auto &shape : other.packing)
        this->addShape(shape->clone());
}

Packing &Packing::operator=(Packing other) {
    std::swap(this->packing, other.packing);
    return *this;
}

Packing::Packing(Packing &&other) noexcept : packing(std::move(other.packing)) {
    other.packing.clear();
}

void Packing::store(std::ostream &out) const {
    for(auto &shape : this->packing)
        shape->store(out);
}

void Packing::restore(std::istream &in) {
    RND rnd(1);
    packing.clear();

    while (!in.eof()) {
        RSAShape *shape = ShapeFactory::createShape(&rnd);
        shape->restore(in);
        this->addShape(shape);
    }
    this->removeShape(this->size() - 1);
}

void Packing::store(const std::string &filename) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open file " + filename + " to store packing");
    this->store(file);
    file.close();
}

void Packing::restore(const std::string &filename) {
    std::ifstream file(filename, std::ios_base::binary);
    if (!file) throw std::runtime_error("Cannot open file " + filename + " to restore packing");
    this->restore(file);
    file.close();
}

void Packing::expandOnPBC(double linearSize, double expandMargin) {
    if (linearSize <= 0.0)
        throw std::runtime_error("linearSize <= 0.0)");
    if (expandMargin <= 0.0 || expandMargin >= 0.5)
        throw std::runtime_error("expandMargin <= 0.0 || expandMargin >= 0.5");

    for (std::size_t i = 0; i < 2; i++) {
        std::size_t oldSize = this->size();
        for (std::size_t j = 0; j < oldSize; j++) {
            auto shape = (*this)[j];
            RSAVector position = shape->getPosition();
            if (position[i] < expandMargin * linearSize)
                expandShapeOnBC(shape, linearSize, i);
            else if (position[i] > (1 - expandMargin) * linearSize)
                expandShapeOnBC(shape, -linearSize, i);
        }
    }
}

/* Helper method. Clones a shape and translates one of its coordinates in a given direction. */
void Packing::expandShapeOnBC(const RSAShape *shape, double translation, size_t translateCoordIdx) {
    RSAShape *shapeClone = shape->clone();
    RSAVector trans;
    trans[translateCoordIdx] = translation;
    shapeClone->translate(trans);
    this->addShape(shapeClone);
}

void Packing::removeShape(std::size_t index) {
    delete (*this)[index];  // free range check
    this->packing.erase(this->packing.begin() + index);
}

const RSAShape *Packing::operator[](std::size_t index) const {
    if (index >= this->size())
        throw std::runtime_error("index >= size");
    return this->packing[index];
}

void Packing::clear() {
    for (auto shape : this->packing)
        delete shape;
    this->packing.clear();
}

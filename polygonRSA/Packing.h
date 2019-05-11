//
// Created by PKua on 17.06.18.
//

#ifndef RSA3D_PACKING_H
#define RSA3D_PACKING_H

#include "shape/Shape.h"
#include <memory>

/**
 * @brief A set of static shapes places on a sufrace.
 *
 * Receives dynamically allocated shapes and take care od deallocation. Works like std::vector with additional specific
 * features.
 */
class Packing {
private:
    using shape_ptr = const RSAShape *;
    using const_iterator = std::vector<shape_ptr>::const_iterator;

    std::vector<shape_ptr> packing;

    void expandShapeOnBC(const RSAShape *shape, double translation, size_t translateCoordIdx);

public:
    Packing() = default;
    virtual ~Packing();

    Packing(const Packing &other);
    Packing(Packing &&other) noexcept;
    Packing &operator=(Packing other);

    /**
     * @breif Adds given dynamically allocated @a shape to a packing and takes care of deallocation later
     * @param shape a Shape to add to a packing
     */
    void addShape(const RSAShape *shape) { this->packing.push_back(shape); }

    /**
     * @brief Access to a shape of given @a index.
     *
     * Internal access to std::vector can throw when index is invalid.
     * @param index index of a shape
     * @return a shape of given @a index
     */
    const RSAShape *operator[](std::size_t index) const;

    /**
     * @brief Removes a shape of given @a index from a packing.
     *
     * Internal access to std::vector can throw when index is invalid.
     * @param index index of a shape
     * @param index
     */
    void removeShape(std::size_t index);

    /**
     * @brief Returns a size of a packing.
     * @return a size of a packing
     */
    std::size_t size() const { return this->packing.size(); }

    /**
     * @brief Removes all shapes from a packing.
     */
    void clear();

    /**
     * @brief Returns const view on internal std::vector which stores shapes.
     * @return const view on internal std::vector which stores shapes
     */
    const std::vector<shape_ptr> &getVector() const { return this->packing; }

    /**
     * @brief Returns the last shape in a packing.
     * @return the last shape in a packing
     */
    const RSAShape *back() const { return this->packing.back(); }

    /**
     * @brief Returns the first shape in a packing.
     * @return the first shape in a packing
     */
    const RSAShape *front() const { return this->packing.front(); }

    /**
     * @brief Stores a packing in a binary form to a given std::ostream.
     * @param out std::ostream to store packing onto
     */
    void store(std::ostream &out) const;

    /**
     * @brief Stores a packing in a binary form to a file of given name.
     *
     * Throws, when a file cannot be accessed.
     * @param filename file name to store a packing to
     */
    void store(const std::string &filename) const;

    /**
     * @brief Restores a packing in a binary form from a given std::istream.
     * @param in std::istream to store packing from
     */
    void restore(std::istream &in);

    /**
     * @brief Restores a packing in a binary form from a file of given name.
     *
     * Throws, when a file cannot be accessed.
     * @param filename file name to restore a packing from
     */
    void restore(const std::string &filename);

    /**
     * @brief Returns const iterator pointing to the beginning.
     * @return const iterator pointing to the beginning
     */
    const_iterator begin() const { return this->packing.begin(); }

    /**
     * @brief Returns const iterator pointing to the end.
     * @return const iterator pointing to the end
     */
    const_iterator end() const { return this->packing.end(); }

    /**
     * @brief Clones and translates shapes distant from the surface boundary not more than @a expandMargin x @a size
     * according to periodic boundary conditions
     * @param packing a packing to expand
     * @param linearSize the size of the packing
     * @param expandMargin the distance (from 0 to 0.5) relative to @a size from the surface boundary from which shapes
     * will be translated
     */
    void expandOnPBC(double linearSize, double expandMargin);
};


#endif //RSA3D_PACKING_H

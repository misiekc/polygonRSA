/*
 * Positioned.h
 *
 *  Created on: 23.03.2017
 *      Author: ciesla
 */

#ifndef POSITIONED_H_
#define POSITIONED_H_

// #include "Vector.h"

#include <array>
#include "Vector.h"

/**
 * @brief An object located in @a SPATIAL_DIMENSION dimensional space.
 *
 * The class abstracts any object which position can be described by @a SPATIAL_DIMENSION numbers. Its position can
 * be accessed or translated by a given vector.
 * @tparam SPATIAL_DIMENSION number of dimensions of a space in which positioned lives
 */
template <unsigned short SPATIAL_DIMENSION>
class Positioned {
	static_assert(SPATIAL_DIMENSION > 0, "SPATIAL_DIMENTION == 0");

private:
    Vector<SPATIAL_DIMENSION> position;

    /* Helper class computing voxel offsets. Usage of std::array gives free range check */
    class Offset {
    private:
        using specific_vertex_offset = std::array<int, SPATIAL_DIMENSION>;

        std::array<specific_vertex_offset, 1 << SPATIAL_DIMENSION> offsets;

    public:
        Offset();

        const specific_vertex_offset &operator[](std::size_t vertex) const;
    };

protected:

	/**
     * @brief Sets new position.
     *
     * Derived shape classes have to override this method when they want to keep track of
     * Positioned position. One would then typically write:
     * \code
     * void Derived::setPosition(const double *position) {
     *     Positioned::setPosition(position);
     *     // invoke getPosition and for example compute vertices
     * }
     * \endcode
     * @param position new position
     */
    virtual void setPosition(const Vector<SPATIAL_DIMENSION> &position);

public:

    /**
     * @brief Offsets of vertices of a voxel.
     *
     * Technically it behaves like int[2 ^ SPATIAL_DIMENSION][SPATIAL_DIMENSION] array. The first index is index of
     * vertex of a voxel, the latter is coordinate index. Coordinates represent a hypercubic voxel of size 1 with the
     * first vertex in the origin.
     */
    static const Offset offset;

	virtual ~Positioned() = default;

	/**
	 * @brief Returns position of a Positioned.
	 * @return position of a Positioned
	 */
	virtual const Vector<SPATIAL_DIMENSION> &getPosition() const;

    /**
     * @brief Translates positioned by a given vector @a v.
     *
     * Default implementation uses setPosition, so it is enough to override only that method to keep track of position.
     * @param v a vector to translate by
     */
	virtual void translate(const Vector<SPATIAL_DIMENSION> &v);
};

using RSAPositioned = Positioned<2>;

#include "Positioned.tpp"

#endif /* POSITIONED_H_ */

//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _VECTOR_H
    #define _VECTOR_H

#include "Matrix.h"
#include <type_traits>


// Class forward declaration for friend operators
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
class Vector;


// Friend operators declaraton for friendship
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
std::ostream & operator<<(std::ostream & _ostr, const Vector<DIM, E> & _v);

template <std::size_t DIM, typename E>
E operator*(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator+(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator-(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator*(const Vector<DIM, E> & _v, E _x);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator*(E _x, const Vector<DIM, E> & _v);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator/(const Vector<DIM, E> & _v, E _x);

template <std::size_t DIM1, std::size_t DIM2, typename E>
Vector<DIM2, E> operator*(const Matrix<DIM2, DIM1, E> & _m, const Vector<DIM1, E> & _v);

template <std::size_t DIM, typename E>
Vector<DIM, E> operator^(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);

template <std::size_t DIM, typename E>
bool operator==(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);

template <std::size_t DIM, typename E>
bool operator!=(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2);



// Class declaration
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
class Vector
{    
    // Inline function defined in class body, rest outside
    //---------------------------------------------------------------------------------------

    // Make each Vector instantiation a friend for inter-size operations (ex. linear trans.)
    //----------------------------------------------------------------------------------------
    template <std::size_t, typename>
    friend class Vector;

    template <std::size_t, std::size_t, typename>
    friend class Matrix;

private:

    Matrix<DIM, 1, E> v;        // Backing matrix
    
    template<bool _COND, typename _T>
    using enabled = typename std::enable_if<_COND, _T>::type;     // Enabled-if type alias
    
    
    // Private access
    //----------------------------------------------------------------------------------------
    E & _get(std::size_t coord)     // Coord reference access without range check
    {
        return this->v._get(coord, 0);
    }

    const E & _get(std::size_t coord) const     // Coord const reference access without range check
    {
        return this->v._get(coord, 0);
    }
    
public:

    // Default, copy and move ctor, assignment operators
    // Nothing special here, let the Matrix class do the magic
    //---------------------------------------------------------------------------------------
    Vector()        // Default
    { }
    
    
    // Other ctors
    //---------------------------------------------------------------------------------------
    Vector(std::array<E, DIM> _coords) : v(Matrix<DIM, 1, E>(_coords))      // std::array initialization
    { }
    
    explicit Vector(const E *_arr) : v(Matrix<DIM, 1, E>(_arr))       // C array initialization
    { }
    
    explicit Vector(E _fill) : v(Matrix<DIM, 1, E>(_fill))      // Fill initialization
    { }
    
    explicit Vector(const Matrix<DIM, 1, E> & _v) : v(_v)       // Matrix<DIM, 1, E> copy-initialization
    { }
    
    explicit Vector(Matrix<DIM, 1, E> && _v) : v(std::move(_v))     // Matrix<DIM, 1, E> move-initializtion
    { }
    
    // Friendship declarations
    //---------------------------------------------------------------------------------------
    friend Vector operator+ <> (const Vector & _v1, const Vector & _v2);    // Addition
    friend Vector operator- <> (const Vector & _v1, const Vector & _v2);    // Subtraction
    friend Vector operator* <> (const Vector & _v, E _x);                   // Scalar multiplication
    friend Vector operator* <> (E _x, const Vector & _v);
    friend Vector operator/ <> (const Vector & _v, E _x);                   // Scalar division
    friend E operator* <> (const Vector & _v1, const Vector & _v2);         // Scalar product
    friend Vector operator^ <>(const Vector & _v1, const Vector & _v2);
    friend bool operator== <> (const Vector & _v1, const Vector & _v2);     // Equality
    friend bool operator!= <> (const Vector & _v1, const Vector & _v2);     // Inequality
    
    template <std::size_t DIM1, std::size_t DIM2, typename _E>
    friend Vector<DIM2, _E> operator* (const Matrix<DIM2, DIM1, _E> & _m, const Vector<DIM1, _E> & _v);     // Linear trans. by matrix
   
    friend std::ostream & operator<< <> (std::ostream & _ostr, const Vector & _v);      // Stream insertion
    
/*
    template <std::size_t _DIM, typename _E>
    friend
    enabled<_DIM == 3, Vector<_DIM, _E>>
    operator^(const Vector<_DIM, _E> & _v1, const Vector<_DIM, _E> & _v2);          // Vector product
*/

    
    // Assignment and othem member operators 
    //---------------------------------------------------------------------------------------
    Vector & operator+=(const Vector & other)   // Addition assignment
    {
        this->v += other.v;
        return *this;
    }

    Vector & operator-=(const Vector & other)   // Subtraction assignment
    {
        this->v -= other.v;
        return *this;
    }
    
    Vector & operator*=(E x)       // Multiplication by scalar assignment
    {
        this->v *= x;
        return *this;
    }
    
    Vector & operator/=(E x)       // Division by scalar assignment
    {
        this->v /= x;
        return *this;
    }

    Vector operator-(void) const        // Unary minus operator
    {
        return Vector(std::move(-this->v));
    }
    
    template <std::size_t _DIM = DIM, typename _E = E>
    enabled<_DIM == 3, Vector> &
    operator^=(const Vector<3, E> & other);      // Vector product assignemnt


    // Norm of vector (for double)
    //---------------------------------------------------------------------------------------
    template <typename _E = E>
    enabled<std::is_same<_E, double>::value, double>
    norm() const
    {
        return sqrt(this->norm2());
    }

    // Norm square of vector (for double)
    //---------------------------------------------------------------------------------------
    template <typename _E = E>
    enabled<std::is_same<_E, double>::value, double>
    norm2() const
    {
        return (*this) * (*this);
    }

    // Returns projection on vector (for double)
    //---------------------------------------------------------------------------------------
    template <typename _E = E>
    enabled<std::is_same<_E, double>::value, Vector>
    projectOn(Vector _axis) const
    {
        return (((*this) * _axis) / (_axis * _axis)) * _axis;
    }

    // Normalized vector
    //---------------------------------------------------------------------------------------
    template <typename _E = E>
    enabled<std::is_same<_E, double>::value, Vector>
    normalized() const
    {
        return *this / this->norm();
    }

    template <std::size_t _DIM = DIM, typename _E = E>
    enabled<_DIM == 2 && std::is_same<_E, double>::value, double>
    angle() {
        return std::atan2(this->v(1, 0), this->v(0, 0));
    }


    // Public access operators
    //---------------------------------------------------------------------------------------
    E & operator[](std::size_t coord)       // Coord reference access
    {
        return this->v(coord, 0);
    }

    const E & operator[](std::size_t coord) const       // Coord const reference access
    {
        return this->v(coord, 0);
    }
    
    // Access operators
    //---------------------------------------------------------------------------------------
    std::size_t getDimension() const
    {
        return DIM;
    }

    void copyToArray(E * array) const {
        v.copyToArray(array);
    }

    std::string toString() const {
        std::ostringstream ostr;
        ostr << *this;
        return ostr.str();
    }
};

// Include implementations
//--------------------------------------------------------------------------------------------
#include "Vector.tpp"


#endif  // _VECTOR_H

//--------------------------------------------------------------------------------------------
// Class representing vectors and operations on them. Template definitions
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------


// Addition
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator+(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return Vector<DIM, E>(std::move(_v1.v + _v2.v));
}

// Subtration
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator-(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return Vector<DIM, E>(std::move(_v1.v - _v2.v));
}


// Multiplication by scalar
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator*(E _x, const Vector<DIM, E> & _v)
{
    return Vector<DIM, E>(std::move(_x * _v.v));
}


// Multiplication by scalar (different operands order)
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator*(const Vector<DIM, E> & _v, E _x)
{
    return Vector<DIM, E>(std::move(_x * _v.v));
}


// Division by scalar
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator/(const Vector<DIM, E> & _v, E _x)
{
    return Vector<DIM, E>(std::move(_v.v / _x));
}


// Scalar product
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
E operator*(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    E prod = 0;
    for (std::size_t i = 0; i < DIM; i++)
        prod += _v1._get(i) * _v2._get(i);
    return prod;
}


// Linear transformation (multiplication of matrix and this)
//--------------------------------------------------------------------------------------------
template <std::size_t DIM1, std::size_t DIM2, typename E>
Vector<DIM2, E> operator*(const Matrix<DIM2, DIM1, E> & _m, const Vector<DIM1, E> & _v)
{
    return Vector<DIM2, E>(std::move(_m * _v.v));
}


// Equality
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
bool operator==(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return (_v1.v == _v2.v);
}


// Inequality
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
bool operator!=(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return (_v1.v != _v2.v);
}

/*
// Cross product
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E>::enabled<DIM == 3, Vector<DIM, E>>
operator^(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return Vector<3, E>({
        _v1._get(1) * _v2._get(2) - _v1._get(2) * _v2._get(1),
        _v1._get(2) * _v2._get(0) - _v1._get(0) * _v2._get(2),
        _v1._get(0) * _v2._get(1) - _v1._get(1) * _v2._get(0) });
}
*/

// Cross product
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
Vector<DIM, E> operator^(const Vector<DIM, E> & _v1, const Vector<DIM, E> & _v2)
{
    return Vector<3, E>({
        _v1._get(1) * _v2._get(2) - _v1._get(2) * _v2._get(1),
        _v1._get(2) * _v2._get(0) - _v1._get(0) * _v2._get(2),
        _v1._get(0) * _v2._get(1) - _v1._get(1) * _v2._get(0) });
}

// Print vector to _ostr output stream
//-----------------------------------------------------------------------------------------
template <std::size_t DIM, typename E>
std::ostream & operator<<(std::ostream & _ostr, const Vector<DIM, E> & _v)
{
    _ostr << "{";
    for (std::size_t i = 0; i < DIM - 1; i++)
        _ostr << _v.v(i, 0) << ", ";
    _ostr << _v.v(DIM - 1, 0) << "}";
    return _ostr;
}

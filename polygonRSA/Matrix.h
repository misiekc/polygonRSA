//-----------------------------------------------------------------------------------------------------------------------------
// Class representing matrices and operations on them
//-----------------------------------------------------------------------------------------------------------------------------
// (C)PKua 2017
//-----------------------------------------------------------------------------------------------------------------------------

#ifndef _MATRIX_H
    #define _MATRIX_H

#include <array>
#include <utility>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>


// Forward declare the Vector class for friendship
//--------------------------------------------------------------------------------------------
template <std::size_t DIM, typename E = double>
class Vector;

// Forward declare the Matrix class for friend operators (ex. == and !=)
//--------------------------------------------------------------------------------------------
template <std::size_t ROWS, std::size_t COLS, typename E = double>
class Matrix;

// Declare operators for friendship
//--------------------------------------------------------------------------------------------

/**
 * @brief Compares two matrices if they are equal.
 * @param _m1 first matrix to compare
 * @param _m2 second matrix to compare
 * @return true, is @a _m1 is equal to @a _m2, false otherwise
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
bool operator==(const Matrix<ROWS, COLS, E> & _m1, const Matrix<ROWS, COLS, E> & _m2);

/**
 * @brief Compares two matrices if they are not equal.
 * @param _m1 first matrix to compare
 * @param _m2 second matrix to compare
 * @return true, is _m1 is not equal to _m2, false otherwise
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
bool operator!=(const Matrix<ROWS, COLS, E> & _m1, const Matrix<ROWS, COLS, E> & _m2);

/**
 * @brief Prints @a _matrix on @a _stream.
 *
 * Output format:
 * @code
 * 1 2 3
 * 4 5 6
 * 7 8 9
 * @endcode
 * @param _stream stream to print matrix on
 * @param _matrix matrix to print
 * @return a reference to @a _stream
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
std::ostream & operator<< (std::ostream & _stream, const Matrix<ROWS, COLS, E> & _matrix);

/**
 * @brief Multiplies @a matrix1 matrix by @a matrix2 an returns the result.
 *
 * Template instantiation will fail if matrices cannot be multiplied.
 * @tparam ROWS1 number of rows in @a matrix1
 * @tparam ROWS_COLS number of columns in @a matrix1 and number of rows in @a matrix2. They must be equal.
 * @tparam COLS2 number of columns in @a matrix2
 * @param matrix1 left operand of multiplication
 * @param matrix2 right operand of multiplication
 * @return result of multiplication
 */
template <std::size_t ROWS1, std::size_t ROWS_COLS, std::size_t COLS2, typename _E>
Matrix<ROWS1, COLS2, _E> operator*(const Matrix<ROWS1, ROWS_COLS, _E> & matrix1, const Matrix<ROWS_COLS, COLS2, _E> & matrix2);

/**
 * @brief Adds @a matrix2 to @a matrix1 and returns result.
 * @param matrix1 left operand matrix
 * @param matrix2 right operand matrix
 * @return @a matrix1 and @a matrix2 sum
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator+(Matrix<ROWS, COLS, E> matrix1, const Matrix<ROWS, COLS, E> & matrix2);

/**
 * @brief Subtracts @a matrix2 from @a matrix1 and returns result.
 * @param matrix1 left operand matrix
 * @param matrix2 right operand matrix
 * @return @a matrix1 and @a matrix2 difference
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator-(Matrix<ROWS, COLS, E> matrix1, const Matrix<ROWS, COLS, E> & matrix2);

/**
 * @brief Multiplies @a matrix by a scalar @a x. Version `matrix * x`.
 * @param matrix matrix to multiply
 * @param x scalar to multiply by
 * @return @a matrix multiplied by a scalar @x
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator*(Matrix<ROWS, COLS, E> matrix, E x);

/**
 * @brief Multiplies @a matrix by a scalar @a x. Version `x * matrix`.
 * @param matrix matrix to multiply
 * @param x scalar to multiply by
 * @return @a matrix multiplied by a scalar @x
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator*(E x, Matrix<ROWS, COLS, E> matrix);

/**
 * @brief divides @a matrix by a scalar @a x.
 * @param matrix matrix to multiply
 * @param x scalar to multiply by
 * @return @a matrix divided by a scalar @x
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator/(Matrix<ROWS, COLS, E> matrix, E x);

/**
 * @brief Unary minus. Returns a matrix consisting of the same elements as `this` but with opposite signs
 * @return a matrix consisting of the same elements as `this` but with opposite signs
 */
template <std::size_t ROWS, std::size_t COLS, typename E>
Matrix<ROWS, COLS, E> operator-(Matrix<ROWS, COLS, E> matrix);


// Matrix class declaration
//--------------------------------------------------------------------------------------------

/**
 * @brief A simple class representing matrices and allowing simple operations on them.
 *
 * A class is best described by some examples of usage:
 * @code
 * Matrix<2, 2> mat1; // Constructs zero-filled 2x2 matrix of doubles
 * mat1(0, 1) = 3; // Sets 3 is a cell from 1st row and second column
 *
 * Matrix<2, 3> mat2{{1, 2, 3, 4, 5, 6}}; // Creates a matrix from std::array, row by row
 *
 * Matrix<2, 3> mat1x2 = mat1 * mat2; // Overloaded operators for basic actions.
 * // Compiler error - cannot multiply 2x3 matrix by 2x2 matrix
 * // Matrix<?, ?> cannot_multiply = mat2 * mat1;
 *
 * // Some more stuff
 * if (Matrix<2, 2>::identity() * mat1 == mat1)
 *     std::cout << "Mathematic works" << std::endl;
 *
 * Matrix<2, 2> rotation = Matrix<2, 2>::rotation(someAngle);
 * Matrix<3, 3> rotation = Matrix<3, 3>::rotation(angleX, angleY, angleZ);
 * @endcode
 * @tparam ROWS number of rows in matrix
 * @tparam COLS number of cols in matrix
 * @tparam E type of element of a matrix. Defaults to `double`
 */
template
<std::size_t ROWS, std::size_t COLS, typename E>
class Matrix
{
private:
    // Make every instantiation friend for inter-size opearations (transpose, minor)
    template <std::size_t, std::size_t, typename>
    friend class Matrix;

    // Allow vector class to access arr
    template <std::size_t, typename>
    friend class Vector;

    E arr[ROWS * COLS];
    
    E & _get(std::size_t row, std::size_t column);
    const E & _get(std::size_t row, std::size_t column) const;

public:
    // Default, copy and move ctor, copy and move assingment operation, dtor
    //----------------------------------------------------------------------------------------

    /// @brief Constructs zero-filled matrix.
    Matrix();
    
    // Other ctors
    //----------------------------------------------------------------------------------------

    /**
     * @brief Constructs a matrix filled with @a _fill value.
     * @param _fill a value to fill matrix
     */
    explicit Matrix(E _fill);

    /**
     * @brief Constructs a matrix from array of pointer to arrays.
     * @param _arr array of pointer to arrays. Successive arrays correspond to matrix rows
     */
    Matrix(const E *const *_arr);

    /**
     * @brief Constructs a matrix from one-dimentional array of matrix elements row by row.
     * @param _arr
     */
    Matrix(const E *_arr);

    /**
     * @brief Constructs a matrix from std::array.
     * @param _arr array to construct from
     */
    Matrix(const std::array<E, ROWS * COLS> & _arr);

    // Operators
    //----------------------------------------------------------------------------------------

    template <std::size_t ROWS1, std::size_t ROWS_COLS, std::size_t COLS2, typename _E>
    friend Matrix<ROWS1, COLS2, _E> operator*(const Matrix<ROWS1, ROWS_COLS, _E> & matrix1, const Matrix<ROWS_COLS, COLS2, _E> & matrix2);

    /**
     * @brief Adds @a other matrix to `this` and assigns result to `this`.
     * @param other matrix to add
     * @return `*this`
     */
    Matrix<ROWS, COLS, E> & operator+=(const Matrix<ROWS, COLS, E> & other);

    /**
     * @brief Subtracts @a other matrix from `this` and assigns result to `this`.
     * @param other matrix to subtract
     * @return `*this`
     */
    Matrix<ROWS, COLS, E> & operator-=(const Matrix<ROWS, COLS, E> & other);

    /**
     * @brief Divides all elements of matrix by @a x and assigns result to `this`.
     * @param x scalar to divide matrix by
     * @return `*this`
     */
    Matrix<ROWS, COLS, E> & operator/=(E x);

    /**
     * @brief Multiplies matrix by @a other, resulting in a matrix of the same size as initial and assigns result to
     * `this`.
     * @param other a matrix to multiply by
     * @return `*this`
     */
    Matrix<ROWS, COLS, E> & operator*=(const Matrix<COLS, COLS, E> & other);

    /**
     * @brief Multiplies all elements of matrix by @a x and assigns result to `this`.
     * @param x scalar to divide matrix by
     * @return `*this`
     */
    Matrix<ROWS, COLS, E> & operator*=(E x);

    /**
     * @brief Gives read/write access to the element from @a row + 1 row and @a column + 1 column.
     *
     * Matrix indices start from 0, like array, not from 1, like in linear algebra.
     * @param row a row index reduces by one
     * @param column a column index reduces by one
     * @return the element from @a row + 1 row and @a column + 1 column
     */
    E & operator()(std::size_t row, std::size_t column);

    /**
     * @brief Gives read-only access to the element from @a row + 1 row and @a column + 1 column.
     *
     * Matrix indices start from 0, like array, not from 1, like in linear algebra.
     * @param row a row index reduces by one
     * @param column a column index reduces by one
     * @return the element from @a row + 1 row and @a column + 1 column
     */
    const E & operator()(std::size_t row, std::size_t column) const;

    /**
     * @brief Returns a row of the matrix as a Vector
     * @param row index of a row (starting from 0) to extract
     * @return a row of the matrix as a Vector
     */
    Vector<COLS, E> row(std::size_t row) const;

    /**
     * @brief Returns a column of the matrix as a Vector
     * @param column index of a column (starting from 0) to extract
     * @return a column of the matrix as a Vector
     */
    Vector<ROWS, E> column(std::size_t column) const;
    
    // Friendship declarations
    //----------------------------------------------------------------------------------------

    friend bool operator== <> (const Matrix & _m1, const Matrix & _m2);
    friend bool operator!= <> (const Matrix & _m1, const Matrix & _m2);
    friend std::ostream & operator<< <> (std::ostream & _stream, const Matrix & _matrix);
    friend Matrix operator+ <> (Matrix matrix1, const Matrix &matrix2);
    friend Matrix operator- <> (Matrix matrix1, const Matrix &matrix2);
    friend Matrix operator* <> (Matrix matrix, E x);
    friend Matrix operator* <> (E x, Matrix matrix);
    friend Matrix operator/ <> (Matrix matrix, E x);
    friend Matrix operator- <> (Matrix matrix);

    // Other operations
    //----------------------------------------------------------------------------------------

    /**
     * @brief Returns number of rows in a matrix.
     * @return number of rows in a matrix
     */
    std::size_t getRows() const;

    /**
     * @brief Returns number of columns in a matrix.
     * @return number of columns in a matrix
     */
    std::size_t getCols() const;

    /**
     * @brief Returns transposition of a matrix.
     * @return transposition of a matrix
     */
    Matrix<COLS, ROWS, E> transpose() const;

    /**
     * @brief Returns string representation of a matrix, as described in operator<<(std::ostream&, const Matrix&).
     * @return string representation of a matrix
     */
    std::string toString() const;
    
    // Identity matrix
    //----------------------------------------------------------------------------------------

    /**
     * @brief Returns identity matrix.
     *
     * Enabled only for n x n matrices, so one shall write
     * @code
     * Matrix<2,2>::identity();
     * Matrix<3,3>::identity();
     * Matrix<10,10>::identity();
     * @endcode
     * @return identity matrix
     */
    template<std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    static typename std::enable_if<_ROWS == _COLS, Matrix<ROWS, COLS, E>>::type 
    identity();

    /**
     * @brief Returns trace of a matrix. Enabled for square matrices.
     *
     * @return trace of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS, E>::type
    tr() const;
    
    // Family of det methods
    //----------------------------------------------------------------------------------------

    /**
     * @brief Returns determinant of a matrix. Version for 1 x 1 matrices.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @return determinant of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 1, E>::type
    det() const;

    /**
     * @brief Returns determinant of a matrix. Version for 2 x 2 matrices.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @return determinant of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
    det() const;

    /**
     * @brief Returns determinant of a matrix. Version for 3 x 3 matrices.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @return determinant of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 3, E>::type
    det() const;

    /**
     * @brief Returns determinant of a matrix. Version for n x n matrices, with n > 3.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @return determinant of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 3), E>::type
    det() const;
    
    // Family of inverse methods
    //----------------------------------------------------------------------------------------

    /**
     * @brief Returns the inverse of a matrix. Version for 1 x 1 matrices.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @throws std::runtime_error if `det()` equals zero
     * @return the inverse of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 1, Matrix<ROWS, COLS, E>>::type
    inverse() const;

    /**
     * @brief Returns the inverse of a matrix. Version for 2 x 2 matrices.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @throws std::runtime_error if det() equals zero
     * @return the inverse of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, Matrix<ROWS, COLS, E>>::type
    inverse() const;

    /**
     * @brief Returns the inverse of a matrix. Version for n x n matrices, with n > 2.
     *
     * Other 3 versions are for matrices of different sizes - compiler will choose appropriate by itself.
     * @throws std::runtime_error if det() equals zero
     * @return the inverse of a matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), Matrix<ROWS, COLS, E>>::type
    inverse() const;
    
    // Family of matrix_minor methods
    //----------------------------------------------------------------------------------------

    /**
     * @brief Returns the minor of a matrix. Version for 2 x 2 matrices.
     *
     * It is calculated by removing @a _row + 1 row and @a column + 1 column and calculating det().
     * @param _row
     * @param _column
     * @return
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && _ROWS == 2, E>::type
    matrix_minor(std::size_t _row, std::size_t _column) const;

    /**
     * @brief Returns the minor of a matrix. Version for n x n matrices, where n > 1.
     *
     * It is calculated by removing @a _row + 1 row and @a column + 1 column and calculating det().
     * @param _row
     * @param _column
     * @return
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS>
    typename std::enable_if<_ROWS == _COLS && (_ROWS > 2), E>::type
    matrix_minor(std::size_t _row, std::size_t _column) const;
    
    // Rotation matrices
    //----------------------------------------------------------------------------------------

    /**
     * @brief Generates 2x2 coutter-clockwise rotation matrix for given angle _a.
     *
     * Enabled only in Matrix<2, 2> so one shall write:
     * \code
     * Matrix<2,2>::rotation(M_PI)
     * \endcode
     * @param _a angle of rotation
     * @return 2x2 rotation matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS, typename _E = E>
    static typename std::enable_if<_ROWS == 2 && _COLS == 2 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
    rotation(double _a);

    /**
     * @brief Generates 3x3 rotation matrix for given angles _ax, _ay, _az.
     *
     * Rotations are counterclockwise and performed around
     * X, Y, Z axis by _ax, _ay, _az respectively in mentioned order. Enabled only in Matrix<3, 3> so one shall write:
     * \code
     * Matrix<3,3>::rotation(0, M_PI/3, M_PI/4)
     * \endcode
     * @param _ax angle of rotation around X axis
     * @param _ay angle of rotation around Y axis
     * @param _az angle of rotation around Z axis
     * @return 3x3 rotation matrix
     */
    template <std::size_t _ROWS = ROWS, std::size_t _COLS = COLS, typename _E = E>
    static typename std::enable_if<_ROWS == 3 && _COLS == 3 && std::is_same<_E, double>::value, Matrix<ROWS, COLS, E>>::type
    rotation(double _ax, double _ay, double _az);

    // Access
    //----------------------------------------------------------------------------------------

    /**
     * @brief Copies matrix elements, row by row, to @a arr array.
     * @param arr a pointer to array to copy matrix to
     */
    inline void copyToArray(E * arr) const;
};



#include "Matrix.tpp"


#endif  // _MATRIX_H

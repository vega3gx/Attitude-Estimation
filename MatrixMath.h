#ifndef MATRIX_MATH_H
#define MATRIX_MATH_H

/**
 * @file
 *
 * @section DESCRIPTION
 *
 * This file implements a basic 3x3 matrix math library. Basic matrix operations are provided along
 * with the matrix inverse function (though that function cannot handle singular matrices).
 *
 * Matrices are defined in row-major order, so that the matrix:
 *   0 1 2
 *   3 4 5
 *   6 7 8
 * is represented by the array `float mat[3][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};`.
 */

/**
 * FP_DELTA defines the tolerance for testing equality for floating-point numbers. Used within
 * MatrixEquals() and is useful for comparison of return values of MatrixDeterminant() and
 * MatrixTrace().
 */
#define FP_DELTA 0.0001

/******************************************************************************
 * Matrix - Matrix Operations
 *****************************************************************************/

/**
 * MatrixMultiply performs a matrix-matrix multiplication operation on two 3x3
 * matrices and returns the result in the third argument.
 */
void MatrixMultiply(double mat1[3][3], double mat2[3][3], double result[3][3]);

/**
 * MatrixAdd performs a matrix addition operation on two 3x3 matrices and 
 * returns the result in the third argument.
 */
void MatrixAdd(double mat1[3][3], double mat2[3][3], double result[3][3]);

/**
 * MatrixEquals checks if the two matrix arguments are equal. The result is
 * 0 if FALSE and 1 if TRUE to follow the standard C conventions of TRUE and
 * FALSE.
 */
int MatrixEquals(double mat1[3][3], double mat2[3][3]);

/******************************************************************************
 * Matrix - Scalar Operations
 *****************************************************************************/

/**
 * MatrixScalarMultiply performs the multiplication of a matrix by a scalar.
 * The result is returned in the third argument.
 */
void MatrixScalarMultiply(double x, double mat[3][3], double result[3][3]);

/**
 * MatrixScalarAdd performs the addition of a matrix by a scalar. The result
 * is returned in the third argument.
 */
void MatrixScalarAdd(double x, double mat[3][3], double result[3][3]);

/******************************************************************************
 * Unary Matrix Operations
 *****************************************************************************/

/**
 * MatrixDeterminant calculates the determinant of a matrix and returns the
 * value as a float.
 */
double MatrixDeterminant(double mat[3][3]);

/**
 * MatrixTrace calculates the trace of a matrix. The result is returned as a
 * float.
 */
double MatrixTrace(double mat[3][3]);

/**
 * MatrixTranspose calculates the transpose of a matrix and returns the
 * result through the second argument
 */
void MatrixTranspose(double mat[3][3], double result[3][3]);

/**
 * MatrixInverse calculates the inverse of a matrix and returns the
 * result through the second argument.
 */
void MatrixInverse(double mat[3][3], double result[3][3]);

void vectorCopy(double source[3], double dest[3]);

void MatrixVectorProduct(double mat[3][3], double vect[3], double prod[3]);

void VectorScale(double vect[3], double scale, double out[3]);

void VectorAdd(double vect1[3], double vect2[3], double sum[3]);

void MatrixCopy(double in[3][3], double out[3][3]);

/**
 * MatrixPrint sends a 3x3 array to standard output with clean formatting.
 * The formatting does not need to look like the expected output given to you
 * in MatricMath.c but each element of the matrix must be separated and have
 * distinguishable position (more than a single line of output).
 */
void MatrixPrint(double mat[3][3]);

void VectorPrint(double vect[3]);

/******************************************************************************
 * Extra Credit Matrix Operations
 *****************************************************************************/

/**
 * Calculate the adjugate matrix of a 3x3 matrix. This is the transpose of the cofactor matrix.
 * @param mat The input 3x3 matrix.
 * @param result The output 3x3 matrix.
 */
void MatrixAdjugate(double mat[3][3], double result[3][3]);

#endif // MATRIX_MATH_H

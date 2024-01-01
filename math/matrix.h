//
// Created by peter on 1/1/2024.
//

#ifndef C_LIB_MATRIX_H
#define C_LIB_MATRIX_H

#include "complex.h"

#define MATRIX_MAX_SIZE 4
#define MATRIX_QR_EIGEN_ITERATION 300

typedef struct {
    int row;
    int col;
    float a[MATRIX_MAX_SIZE][MATRIX_MAX_SIZE];
} Matrix;

struct MatrixControl {
    //////////////////////////BASIC//////////////////////////
    /**
     * This will give matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param a element of matrix in 2D array of  number
     * @return matrix
     */
    Matrix (*get)(int row, int col, float a[row][col]);

    /**
     * This will give matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param c  number for all the elements in matrix
     * @return matrix x+iy
     */
    Matrix (*getConstantMatrix)(int row, int col, float c);

    //////////////////////////////SPECIAL MATRIX///////////////////////////
    /**
     * Gives identity matrix
     * @param size size of matrix
     * @return identity matrix
     */
    Matrix (*identity)(int size);

    /**
     * Gives null matrix
     * @param row row of null matrix
     * @param col column of null matrix
     * @return null matrix
     */
    Matrix (*null)(int row, int col);

    /////////////////////////////ARITHMETIC OPERATIONS/////////////////////
    /**
     * Add two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 + m2
     */
    Matrix (*add)(Matrix m1, Matrix m2);

    /**
     * Gives differences of two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 - m2
     */
    Matrix (*sub)(Matrix m1, Matrix m2);

    /**
     * Gives product of two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 x m2
     */
    Matrix (*dot)(Matrix m1, Matrix m2);

    /**
     * Gives scaled matrix
     * @param m matrix
     * @param s scale factor
     * @return s x m
     */
    Matrix (*scale)(Matrix m, float s);

    /////////////////////////MATRIX OPERATIONS/////////////////////////
    /**
     * Gives minor matrix
     * @param m matrix whose minor is to be given
     * @param i row index
     * @param j column index
     * @return minor of matrix m at (i,j)
     */
    Matrix (*minor)(Matrix m, int i, int j);

    /**
     * Gives determinant of matrix
     * @param m given matrix
     * @return determinant of matrix m
     */
    float (*det)(Matrix m);

    /**
     * Gives cofactor of matrix
     * @param m matrix
     * @param i row index
     * @param j column index
     * @return cofactor of m at (i,j)
     */
    float (*cofactor)(Matrix m, int i, int j);

    /**
     * Gives cofactor matrix
     * @param m given matrix
     * @return cofactor matrix of m
     */
    Matrix (*cofactorMatrix)(Matrix m);

    /**
     * Gives transpose of matrix
     * @param m given matrix
     * @return m' (transpose matrix is not conjugate matrix!!!)
     */
    Matrix (*transpose)(Matrix m);

    /**
     * Gives adjoint matrix
     * @param m given matrix
     * @return transpose of cofactor of m
     */
    Matrix (*adjoint)(Matrix m);

    /**
     * Gives inverse of matrix
     * @param m given matrix
     * @return adjoint(m)/det(m)
     */
    Matrix (*inverse)(Matrix m);

    /**
     * Gives normalize matrix
     * @param m given matrix
     * @return normalize each column or vector in matrix m
     */
    Matrix (*normalise)(Matrix m);

    /**
     * Gives rank of matrix
     * @param m given matrix
     * @return rank of m
     */
    int (*rank)(Matrix m);

    /////////////////////////DECOMPOSITION//////////////////////////
    /**
     * Gives Given Rotation's matrix. this makes (i2,j) element 0 after operation on given matrix
     * @param m given matrix
     * @param i1 first row
     * @param i2 second row
     * @return  given rotation about row=col=i1 with span i2-i1
     */
    Matrix (*givenRotationMatrix)(Matrix m, int i1, int i2);

    /**
     * This QR decomposition using Given rotation matrix (m = Q x R)
     * @param m given matrix
     * @param Q Orthogonal matrix
     * @param R Right triangular matrix
     */
    void (*getQRWithGiven)(Matrix m, Matrix *Q, Matrix *R);

    /**
     * Gives QR with Gram - Schmidt method (m = Q x R)
     * @param m given matrix
     * @param Q Orthogonal matrix
     * @param R Right triangular matrix
     */
    void (*getQRWithGram)(Matrix m, Matrix *Q, Matrix *R);

    /**
     * Check if given matrix is identity
     * @param m given matrix
     * @return returns 1 if m is identity
     */
    int (*isIdentityMatrix)(Matrix m);

    /**
     * Check if given matrix is diagonal matrix
     * @param m given matrix
     * @return returns 1 if m is diagonal matrix
     */
    int (*isDiagonalMatrix)(Matrix m);

    /**
     * Check if given matrix is null matrix
     * @param m given matrix
     * @return returns 1 if m is null matrix
     */
    int (*isNullMatrix)(Matrix m);

    /**
     * Gives singular value decomposition of given matrix (m = U * S x conjugate(V))
     * @param m given matrix of size (r x c)
     * @param U Unitary matrix of size (r x r)
     * @param S Rectangular diagonal matrix of size (r x c)
     * @param V Unitary matrix of size (c x c)
     */
    void (*getSVDDecompose)(Matrix m, Matrix *U, Matrix *S, Matrix *V);

    ///////////////////////////OPERATIONS///////////////////////////////////////////////////
    /**
     * Returns matrix applying function f to each diagonal element
     * @param f function that take  number and return  number
     * @param m give matrix
     * @return modified matrix
     */
    Matrix (*diagonalOperation)(float (*f)(float x), Matrix m);

    /**
     * Returns matrix applying function f to each element
     * @param f function that takes  number and return  number
     * @param m give matrix
     * @return modified matrix
     */
    Matrix (*elementWiseOperation)(float (*f)(float x), Matrix m);

    /**
     * This returns matrix after applying analytical function (should be taylor expansion!!!)
     * @param f analytic function that takes  number and return  number
     * @param m give matrix
     * @return modified matrix
     */
    Matrix (*analyticFunction)(Complex (*f)(Complex x), Matrix m);

    /**
     * This returns matrix after applying analytical function (should be taylor expansion!!!)
     * @param f analytic function that takes two  number and return  number
     * @param m given matrix
     * @param c given  number
     * @return modified matrix
     */
    Matrix (*analyticFunction2)(Complex (*f)(Complex x, Complex c), Matrix m, float c);

    //////////////////////////EXPONENTIAL, LOGARITHMIC AND POWER/////////////////////////////
    /**
     * Gives exponential of matrix
     * @param m given matrix
     * @return exp(m)
     */
    Matrix (*exp)(Matrix m);

    /**
     * Gives natural log of matrix
     * @param m given matrix
     * @return ln(m)
     */
    Matrix (*ln)(Matrix m);

    /**
     * Gives pow of matrix raise to  number
     * @param m given matrix
     * @param c given  number
     * @return m^c
     */
    Matrix (*pow)(Matrix m, float c);

    ////////////////////////PRINTING////////////////////////////////////
    /**
     * Print the matrix in form of x + iy
     * @param m given matrix
     */
    void (*print)(Matrix m);

};

__attribute__((unused)) extern struct MatrixControl StaticMatrix;
#endif //C_LIB_MATRIX_H

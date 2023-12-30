//
// Created by peter on 12/29/2023.
//

#ifndef C_LIB_COMPLEX_MATRIX_H
#define C_LIB_COMPLEX_MATRIX_H
#include "complex.h"

#define COMPLEX_MATRIX_MAX_SIZE 4
#define COMPLEX_MATRIX_QR_EIGEN_ITERATION 200

typedef struct {
    int row;
    int col;
    Complex a[COMPLEX_MATRIX_MAX_SIZE][COMPLEX_MATRIX_MAX_SIZE];
} ComplexMatrix;

struct ComplexMatrixControl{
    //////////////////////////BASIC//////////////////////////
    /**
     * This will give complex matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param a element of matrix in 2D array of complex number
     * @return complex matrix
     */
    ComplexMatrix (*get)(int row, int col, Complex a[row][col]);

    /**
     * This will give complex matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param a element of matrix in 2D array of real number
     * @return complex matrix
     */
    ComplexMatrix (*getFromReal)(int row, int col, float a[row][col]);

    /**
     * This will give complex matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param x element of matrix in 2D array of real number
     * @param y element of matrix in 2D array of image number
     * @return complex matrix x+iy
     */
    ComplexMatrix (*getXY)(int row, int col, float x[row][col],float y[row][col]);

    /**
     * This will give complex matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param c complex number for all the elements in matrix
     * @return complex matrix x+iy
     */
    ComplexMatrix (*getConstantMatrix)(int row, int col, Complex c);

    /**
     * This will give complex matrix
     * @param row row of matrix
     * @param col column of matrix
     * @param x real value for all the elements in matrix
     * @return complex matrix x+iy
     */
    ComplexMatrix (*getConstantMatrixFromReal)(int row, int col, float x);

    //////////////////////////////SPECIAL MATRIX///////////////////////////
    /**
     * Gives identity matrix
     * @param size size of matrix
     * @return identity matrix
     */
    ComplexMatrix (*identity)(int size);

    /**
     * Gives null matrix
     * @param row row of null matrix
     * @param col column of null matrix
     * @return null matrix
     */
    ComplexMatrix (*null)(int row, int col);

    /////////////////////////////ARITHMETIC OPERATIONS/////////////////////
    /**
     * Add two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 + m2
     */
    ComplexMatrix (*add)(ComplexMatrix m1, ComplexMatrix m2);

    /**
     * Gives differences of two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 - m2
     */
    ComplexMatrix (*sub)(ComplexMatrix m1, ComplexMatrix m2);

    /**
     * Gives product of two matrices
     * @param m1 first matrix
     * @param m2 second matrix
     * @return m1 x m2
     */
    ComplexMatrix (*dot)(ComplexMatrix m1, ComplexMatrix m2);

    /**
     * Gives scaled matrix
     * @param m matrix
     * @param s scale factor
     * @return s x m
     */
     ComplexMatrix (*scale)(ComplexMatrix m, Complex s);

    /////////////////////////MATRIX OPERATIONS/////////////////////////
    /**
     * Gives minor matrix
     * @param m matrix whose minor is to be given
     * @param i row index
     * @param j column index
     * @return minor of matrix m at (i,j)
     */
    ComplexMatrix (*minor)(ComplexMatrix m, int i, int j);

    /**
     * Gives determinant of matrix
     * @param m given matrix
     * @return determinant of matrix m
     */
    Complex (*det)(ComplexMatrix m);

    /**
     * Gives cofactor of matrix
     * @param m matrix
     * @param i row index
     * @param j column index
     * @return cofactor of m at (i,j)
     */
    Complex (*cofactor)(ComplexMatrix m, int i, int j);

    /**
     * Gives cofactor matrix
     * @param m given matrix
     * @return cofactor matrix of m
     */
    ComplexMatrix (*cofactorMatrix)(ComplexMatrix m);

    /**
     * Gives transpose of matrix
     * @param m given matrix
     * @return m' (transpose matrix is not conjugate matrix!!!)
     */
    ComplexMatrix (*transpose)(ComplexMatrix m);

    /**
     * Gives adjoint matrix
     * @param m given matrix
     * @return transpose of cofactor of m
     */
    ComplexMatrix (*adjoint)(ComplexMatrix m);

    /**
     * Gives inverse of matrix
     * @param m given matrix
     * @return adjoint(m)/det(m)
     */
    ComplexMatrix (*inverse)(ComplexMatrix m);

    /**
     * Gives normalize matrix
     * @param m given matrix
     * @return normalize each column or vector in matrix m
     */
    ComplexMatrix (*normalise)(ComplexMatrix m);

    /**
     * Gives conjugate matrix
     * @param m given matrix
     * @return conjugate each element and transpose the resultant matrix
     */
    ComplexMatrix (*conjugate)(ComplexMatrix m);

    /**
     * Gives rank of matrix
     * @param m given matrix
     * @return rank of m
     */
    int (*rank)(ComplexMatrix m);

    /////////////////////////SPECIAL OPERATIONS//////////////////////////
    /**
     * Gives Given Rotation's matrix. this makes (i2,j) element 0 after operation on given matrix
     * @param m given matrix
     * @param i1 first row
     * @param i2 second row
     * @return  given rotation about row=col=i1 with span i2-i1
     */
    ComplexMatrix (*givenRotationMatrix)(ComplexMatrix m, int i1, int i2);

    /**
     * This QR decomposition using Given rotation matrix (m = Q x R)
     * @param m given matrix
     * @param Q Orthogonal matrix
     * @param R Right triangular matrix
     */
    void (*getQRWithGiven)(ComplexMatrix m, ComplexMatrix *Q, ComplexMatrix *R);

    /**
     * Gives QR with Gram - Schmidt method (m = Q x R)
     * @param m given matrix
     * @param Q Orthogonal matrix
     * @param R Right triangular matrix
     */
    void (*getQRWithGram)(ComplexMatrix m, ComplexMatrix *Q, ComplexMatrix *R);

    /**
     * Check if given matrix is identity
     * @param m given matrix
     * @return returns 1 if m is identity
     */
    int (*isIdentityMatrix)(ComplexMatrix m);

    /**
     * Check if given matrix is diagonal matrix
     * @param m given matrix
     * @return returns 1 if m is diagonal matrix
     */
    int (*isDiagonalMatrix)(ComplexMatrix m);

    /**
     * Check if given matrix is null matrix
     * @param m given matrix
     * @return returns 1 if m is null matrix
     */
    int (*isNullMatrix)(ComplexMatrix m);

    /**
     * Gives eigen decomposition of given matrix (m = V x D)
     * @param m given matrix
     * @param V eigen vector matrix
     * @param D eigen value matrix
     */
     void (*getEigenDecompose)(ComplexMatrix m, ComplexMatrix *V, ComplexMatrix *D);

    /**
     * Gives eigen decomposition of given matrix using QR method(m = V x D)
     * @param m given matrix
     * @param V eigen vector matrix
     * @param D eigen value matrix
     */
    void (*getEigenDecomposeQR)(ComplexMatrix m, ComplexMatrix *V, ComplexMatrix *D);

    /**
     * Gives singular value decomposition of given matrix (m = U * S x conjugate(V))
     * @param m given matrix of size (r x c)
     * @param U Unitary matrix of size (r x r)
     * @param S Rectangular diagonal matrix of size (r x c)
     * @param V Unitary matrix of size (c x c)
     */
    void (*getSVDDecompose)(ComplexMatrix m, ComplexMatrix *U, ComplexMatrix *S, ComplexMatrix *V);

    /**
     * Gives singular value decomposition of given matrix using QR method (m = U * S x conjugate(V))
     * @param m given matrix of size (r x c)
     * @param U Unitary matrix of size (r x r)
     * @param S Rectangular diagonal matrix of size (r x c)
     * @param V Unitary matrix of size (c x c)
     */
    void (*getSVDDecomposeQR)(ComplexMatrix m, ComplexMatrix *U, ComplexMatrix *S, ComplexMatrix *V);

    ////////////////////////PRINTING////////////////////////////////////
    /**
     * Print the matrix in form of x + iy
     * @param m given matrix
     */
    void (*print)(ComplexMatrix m);

    /**
     * Print the matrix in form of r < theta
     * @param m given matrix
     */
    void (*printA)(ComplexMatrix m);
};

__attribute__((unused)) extern struct ComplexMatrixControl StaticComplexMatrix;

#endif //C_LIB_COMPLEX_MATRIX_H

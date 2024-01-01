//
// Created by peter on 1/1/2024.
//

#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include "polynomial.h"
#include "linear.h"

#define min(a, b) (a<=b?a:b)
#define cm StaticComplexMatrix
#define c StaticComplex

//////////////////////////BASIC//////////////////////////////////////////
/**
 * This will give matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param a element of matrix in 2D array of number
 * @return matrix
 */
static Matrix get(int row, int col, float a[row][col]) {
    Matrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = a[i][j];
    return matrix;
}

/**
 * This will give matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param c1  number for all the elements in matrix
 * @return matrix x+iy
 */
static Matrix getConstantMatrix(int row, int col, float c1) {
    Matrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c1;
    return matrix;
}

//////////////////////////////SPECIAL MATRIX///////////////////////////
/**
 * Gives identity matrix
 * @param size size of matrix
 * @return identity matrix
 */
static Matrix identity(int size) {
    Matrix matrix = {.row=size, .col=size};
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix.a[i][j] = i == j ? 1 : 0;
    return matrix;
}

/**
 * Gives null matrix
 * @param row row of null matrix
 * @param col column of null matrix
 * @return null matrix
 */
static Matrix null(int row, int col) {
    Matrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = 0;
    return matrix;
}

/////////////////////////////ARITHMETIC OPERATIONS/////////////////////
/**
 * Add two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 + m2
 */
static Matrix add(Matrix m1, Matrix m2) {
    Matrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            m.a[i][j] = m1.a[i][j]+ m2.a[i][j];
    return m;
}

/**
 * Gives differences of two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 - m2
 */
static Matrix sub(Matrix m1, Matrix m2) {
    Matrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            m.a[i][j] = m1.a[i][j]- m2.a[i][j];
    return m;
}

/**
 * Gives product of two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 x m2
 */
static Matrix dot(Matrix m1, Matrix m2) {
    Matrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            for (int k = 0; k < m1.col; k++)
                m.a[i][j] =m.a[i][j]+ m1.a[i][k]* m2.a[k][j];
    return m;
}

/**
 * Gives scaled matrix
 * @param m matrix
 * @param s scale factor
 * @return s x m
 */
static Matrix scale(Matrix m, float s) {
    Matrix matrix = m;
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            matrix.a[i][j] = matrix.a[i][j]* s;
    return matrix;
}

//////////////////////////////MATRIX OPERATIONS/////////////////////////
/**
 * Gives minor matrix
 * @param m matrix whose minor is to be given
 * @param i row index
 * @param j column index
 * @return minor of matrix m at (i,j)
 */
static Matrix minor(Matrix m, int i, int j) {
    Matrix matrix = null(m.row - 1, m.col - 1);
    for (int k = 0, p = 0; k < m.row; k++) {
        if (k == i)
            continue;
        for (int l = 0, q = 0; l < m.col; l++) {
            if (l == j)
                continue;
            matrix.a[p][q] = m.a[k][l];
            q++;
        }
        p++;
    }
    return matrix;
}

/**
 * Gives determinant of matrix
 * @param m given matrix
 * @return determinant of matrix m
 */
static float det(Matrix m) {
    if (m.row == 1)
        return m.a[0][0];
    float sum = 0;
    for (int j = 0; j < m.col; j++)
        sum = sum+(m.a[0][j]* det(minor(m, 0, j))*powf(-1.0f,(float) j));
    return sum;
}

/**
 * Gives cofactor of matrix
 * @param m matrix
 * @param i row index
 * @param j column index
 * @return cofactor of m at (i,j)
 */
static float cofactor(Matrix m, int i, int j) {
    return (powf(-1.0f, (float) (i + j))*det(minor(m, i, j)));
}

/**
 * Gives cofactor matrix
 * @param m given matrix
 * @return cofactor matrix of m
 */
static Matrix cofactorMatrix(Matrix m) {
    Matrix matrix = null(m.row, m.col);
    for (int i = 0; i < matrix.row; i++)
        for (int j = 0; j < matrix.col; j++)
            matrix.a[i][j] = cofactor(m, i, j);
    return matrix;
}

/**
 * Gives transpose of matrix
 * @param m given matrix
 * @return m' (transpose matrix is not conjugate matrix!!!)
 */
static Matrix transpose(Matrix m) {
    Matrix matrix = null(m.col, m.row);
    for (int i = 0; i < matrix.row; i++)
        for (int j = 0; j < matrix.col; j++)
            matrix.a[i][j] = m.a[j][i];
    return matrix;
}

/**
 * Gives adjoint matrix
 * @param m given matrix
 * @return transpose of cofactor of m
 */
static Matrix adjoint(Matrix m) {
    return transpose(cofactorMatrix(m));
}

/**
 * Gives inverse of matrix
 * @param m given matrix
 * @return adjoint(m)/det(m)
 */
static Matrix inverse(Matrix m) {
    if (m.row == 1 && m.col == 1) {
        Matrix matrix = {.row=1, .col=1, 1.0f/m.a[0][0]};
        return matrix;
    }
    return scale(adjoint(m), 1.0f/det(m));
}

/**
 * Gives normalize matrix
 * @param m given matrix
 * @return normalize each column or vector in matrix m
 */
static Matrix normalise(Matrix m) {
    Matrix matrix = m;
    for (int j = 0; j < matrix.col; j++) {
        float sum = 0;
        for (int i = 0; i < matrix.row; i++)
            sum = sum+ matrix.a[i][j]* matrix.a[i][j];
        sum = powf(sum,0.5f);
        for (int i = 0; i < matrix.row; i++)
            matrix.a[i][j] = matrix.a[i][j]/sum;
    }
    return matrix;
}

/**
 * Gives rank of matrix
 * @param m given matrix
 * @return rank of m
 */
static int rank(Matrix m) {
    Matrix A;
    if (m.row <= m.col)
        A = dot(m, transpose(m));
    else
        A = dot(transpose(m), m);
    A = normalise(A);
    for (int rank = A.row; rank > 0; rank--) {
        if (det(A) < 0.001f) {
            Matrix temp = identity(A.row - 1);
            for (int i = 0; i < temp.row; i++)
                for (int j = 0; j < temp.col; ++j)
                    temp.a[i][j] = A.a[i][j];
            A = temp;
            continue;
        }
        return rank;
    }
    return 0;
}


/////////////////////////SPECIAL OPERATIONS//////////////////////////
/**
 * Gives Given Rotation's matrix. this makes (i2,j) element 0 after operation on given matrix
 * @param m given matrix
 * @param i1 first row
 * @param i2 second row
 * @return  given rotation about row=col=i1 with span i2-i1
 */
static Matrix givenRotationMatrix(Matrix m,int i1, int i2) {
    int j = i1;
    Matrix matrix = identity(m.row);
    float r = powf((m.a[i1][j]* m.a[i1][j])+(m.a[i2][j]* m.a[i2][j]),0.5f);
    int diff = i2 - i1;
    float cos = m.a[i1][j]/r;
    float sine = m.a[i2][j]/r;
    matrix.a[i1][j] = cos;
    matrix.a[i1][j + diff] = sine;
    matrix.a[i2][j] = -sine;
    matrix.a[i2][j + diff] = cos;
    return matrix;
}

/**
 * This QR decomposition using Given rotation matrix (m = Q x R)
 * @param m given matrix
 * @param Q Orthogonal matrix
 * @param R Right triangular matrix
 */
static void getQRWithGiven(Matrix m, Matrix *Q, Matrix *R) {
    *R = m;
    *Q = identity(m.row);
    for (int j = 0; j < m.col; j++) {
        for (int i = 0; i < m.col; i++) {
            if (i <= j)
                continue;
            Matrix G = givenRotationMatrix(*R, j, i);
            *R = dot(G, *R);
            *Q = dot(G, *Q);
        }
    }
    *Q = transpose(*Q);
}

/**
 * Gives QR with Gram - Schmidt method (m = Q x R)
 * @param m given matrix
 * @param Q Orthogonal matrix
 * @param R Right triangular matrix
 */
static void getQRWithGram(Matrix m, Matrix *Q, Matrix *R) {
    //TODO not working with  number
    *Q = null(m.row,m.row);
    *R = null(m.row,m.row);

    Matrix mt = transpose(m);
    Matrix a_t[mt.col];
    Matrix u_n[mt.col];
    for (int i = 0; i < mt.row; i++) {
        a_t[i] = null(1,mt.col);
        for (int j = 0; j < mt.col; ++j)
            a_t[i].a[0][j]=mt.a[i][j];
        u_n[i] = null(mt.row, 1);
    }

    for (int k = 0; k < mt.col; k++) {
        Matrix proj = scale(u_n[0],dot(a_t[k],u_n[0]).a[0][0]);
        u_n[k]=sub(transpose(a_t[k]),proj);
        for (int t = 1; t < k-1; ++t) {
            proj = scale(u_n[t],dot(u_n[k],u_n[t]).a[0][0]);
            u_n[k] = sub(transpose(u_n[k]),proj);
        }
        u_n[k] = normalise(u_n[k]);
    }

    for (int i = 0; i < R->row; i++) {
        for (int j = i; j < R->col; j++)
            R->a[i][j] = dot(a_t[j],u_n[i]).a[0][0];
        for (int j = 0; j < R->col; ++j)
            Q->a[i][j] = u_n[j].a[i][0];
    }
}

/**
 * Check if given matrix is identity
 * @param m given matrix
 * @return returns 1 if m is identity
 */
static int isIdentityMatrix(Matrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if ((i == j) && ((m.a[i][j]-1) > COMPLEX_EPSILON)) {
                return 0;
            } else if (m.a[i][j] > COMPLEX_EPSILON)
                return 0;
    return 1;
}

/**
 * Check if given matrix is diagonal matrix
 * @param m given matrix
 * @return returns 1 if m is diagonal matrix
 */
static int isDiagonalMatrix(Matrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if ((i != j) && (m.a[i][j] > COMPLEX_EPSILON))
                return 0;
    return 1;
}

/**
 * Check if given matrix is null matrix
 * @param m given matrix
 * @return returns 1 if m is null matrix
 */
static int isNullMatrix(Matrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if (m.a[i][j] > COMPLEX_EPSILON)
                return 0;
    return 1;
}


static float eigenCharacteristicsFunction(float x, void **args) {
    Matrix *m = (Matrix *) args[0];
    return det(add(*m, scale(identity(m->row), x)));
}

/**
 * Gives singular value decomposition of given matrix (m = U * S x conjugate(V))
 * @param m given matrix of size (r x c)
 * @param U Unitary matrix of size (r x r)
 * @param S Rectangular diagonal matrix of size (r x c)
 * @param V Unitary matrix of size (c x c)
 */
static void getSVDDecompose(Matrix m, Matrix *U, Matrix *S, Matrix *V) {
    Matrix A = m;
    if (A.row >= A.col) {
        Matrix AAt = dot(A, transpose(A));

        ////////////LEFT SINGULAR MATRIX//////////////////
        Matrix D = AAt;
        Matrix Q,R;
        *U = identity(D.row);
        for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; ++i) {
            getQRWithGiven(D,&Q,&R);
            D = dot(R,Q);
            *U = dot(*U,Q);
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = powf(D.a[i][i], 0.5f);

        ////////////RIGHT SINGULAR MATRIX//////////////////
        Matrix S_inv = transpose(*S);
        for (int i = 0; i < min(S->row, S->col); i++)
            S_inv.a[i][i] = 1/S_inv.a[i][i];
        *V = dot(S_inv, dot(transpose(*U), A));
        *V = transpose(*V);
    } else {
        Matrix AtA = dot(transpose(A), A);
        ////////////RIGHT SINGULAR MATRIX//////////////////
        Matrix D = AtA;
        Matrix Q,R;
        *V = identity(D.row);
        for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; ++i) {
            getQRWithGiven(D,&Q,&R);
            D = dot(R,Q);
            *V = dot(*V,Q);
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = powf(D.a[i][i],0.5f);

        ////////////LEFT SINGULAR MATRIX//////////////////
        Matrix S_inv = transpose(*S);
        for (int i = 0; i < S_inv.col; i++)
            S_inv.a[i][i] = 1.0f/S_inv.a[i][i];
        *U = dot(A, dot(*V, S_inv));
    }
}

///////////////////////////OPERATIONS///////////////////////////////////////////////////
/**
 * Returns matrix applying function f to each diagonal element
 * @param f function that take  number and return  number
 * @param m give matrix
 * @return modified matrix
 */
static Matrix diagonalOperation(float (*f)(float x), Matrix m){
    for (int i = 0; i < min(m.row,m.col); ++i)
        m.a[i][i] = f(m.a[i][i]);
    return m;
}

/**
 * Returns matrix applying function f to each element
 * @param f function that takes  number and return  number
 * @param m give matrix
 * @return modified matrix
 */
static Matrix elementWiseOperation(float (*f)(float x), Matrix m){
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            m.a[i][j] = f(m.a[i][j]);
    return m;
}

/**
 * This returns matrix after applying analytical function (should be taylor expansion!!!)
 * @param f analytic function that takes  number and return  number
 * @param m give matrix
 * @return modified matrix
 */
static Matrix analyticFunction(Complex (*f)(Complex x),Matrix m){
    return cm.toMatrix(cm.analyticFunction(f,cm.getFromMatrix(m)));
}

/**
 * This returns matrix after applying analytical function (should be taylor expansion!!!)
 * @param f analytic function that takes two  number and return  number
 * @param m given matrix
 * @param c1 given  number
 * @return modified matrix
 */
static Matrix analyticFunction2(Complex (*f)(Complex x,Complex c1),Matrix  m,float c1){
    return cm.toMatrix(cm.analyticFunction2(f,cm.getFromMatrix(m),c.getFromReal(c1)));
}

////////////////////////////EXPONENTIAL, LOGARITHMIC AND POWER/////////////////////////////
/**
 * Gives exponential of matrix
 * @param m given matrix
 * @return exp(m)
 */
static Matrix expm(Matrix m){
    return cm.toMatrix(cm.exp(cm.getFromMatrix(m)));
}

/**
 * Gives natural log of matrix
 * @param m given matrix
 * @return ln(m)
 */
static Matrix lnm(Matrix m){
    return cm.toMatrix(cm.ln(cm.getFromMatrix(m)));
}

/**
 * Gives pow of matrix raise to  number
 * @param m given matrix
 * @param c given  number
 * @return m^c
 */
static Matrix powm(Matrix m,float c1){
    return cm.toMatrix(cm.pow(cm.getFromMatrix(m),c.getFromReal(c1)));
}

//////////////////////////PRINTING////////////////////////////////////
/**
 * Print the matrix in form of x + iy
 * @param m given matrix
 */
static void print(Matrix m) {
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j)
            printf("%.3f    ",m.a[i][j]);
        printf("\n");
    }
}


__attribute__((unused)) struct MatrixControl StaticMatrix = {
        .get = get,
        .getConstantMatrix = getConstantMatrix,
        .identity = identity,
        .null = null,
        .add = add,
        .sub = sub,
        .dot = dot,
        .scale = scale,
        .minor = minor,
        .det = det,
        .cofactor = cofactor,
        .cofactorMatrix = cofactorMatrix,
        .transpose = transpose,
        .adjoint = adjoint,
        .inverse = inverse,
        .normalise = normalise,
        .rank = rank,
        .givenRotationMatrix = givenRotationMatrix,
        .getQRWithGiven = getQRWithGiven,
        .getQRWithGram = getQRWithGram,
        .isIdentityMatrix = isIdentityMatrix,
        .isDiagonalMatrix = isDiagonalMatrix,
        .isNullMatrix = isNullMatrix,
        .getSVDDecompose = getSVDDecompose,
        .diagonalOperation = diagonalOperation,
        .elementWiseOperation = elementWiseOperation,
        .analyticFunction = analyticFunction,
        .analyticFunction2 = analyticFunction2,
        .exp = expm,
        .ln = lnm,
        .pow = powm,
        .print = print,
};
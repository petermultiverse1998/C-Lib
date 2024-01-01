//
// Created by peter on 12/29/2023.
//

#include <stdio.h>
#include "complex_matrix.h"
#include "polynomial.h"
#include "linear.h"

#define min(a, b) (a<=b?a:b)
#define c StaticComplex

//////////////////////////BASIC//////////////////////////////////////////
/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param a element of matrix in 2D array of complex number
 * @return complex matrix
 */
static ComplexMatrix get(int row, int col, Complex a[row][col]) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = a[i][j];
    return matrix;
}

/**
 * This will give complex matrix
 * @param m real matrix
 * @return complex matrix
 */
static ComplexMatrix getFromMatrix(Matrix m) {
    ComplexMatrix matrix = {.row=m.row, .col=m.col};
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            matrix.a[i][j] = c.getFromReal(m.a[i][j]);
    return matrix;
}

/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param a element of matrix in 2D array of real number
 * @return complex matrix
 */
static ComplexMatrix getFromReal(int row, int col, float a[row][col]) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c.getFromReal(a[i][j]);
    return matrix;
}

/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param x element of matrix in 2D array of real number
 * @param y element of matrix in 2D array of image number
 * @return complex matrix x+iy
 */
static ComplexMatrix getXY(int row, int col, float x[row][col], float y[row][col]) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c.get(x[i][j], y[i][j]);
    return matrix;
}

/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param c1 complex number for all the elements in matrix
 * @return complex matrix x+iy
 */
static ComplexMatrix getConstantMatrix(int row, int col, Complex c1) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c1;
    return matrix;
}

/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param x real value for all the elements in matrix
 * @return complex matrix x+iy
 */
static ComplexMatrix getConstantMatrixFromReal(int row, int col, float x) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c.getFromReal(x);
    return matrix;
}

//////////////////////////////SPECIAL MATRIX///////////////////////////
/**
 * Gives identity matrix
 * @param size size of matrix
 * @return identity matrix
 */
static ComplexMatrix identity(int size) {
    ComplexMatrix matrix = {.row=size, .col=size};
    Complex c1 = {1, 0};
    Complex c0 = {0, 0};
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix.a[i][j] = i == j ? c1 : c0;
    return matrix;
}

/**
 * Gives null matrix
 * @param row row of null matrix
 * @param col column of null matrix
 * @return null matrix
 */
static ComplexMatrix null(int row, int col) {
    ComplexMatrix matrix = {.row=row, .col=col};
    Complex c0 = {0, 0};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c0;
    return matrix;
}

/////////////////////////////ARITHMETIC OPERATIONS/////////////////////
/**
 * Add two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 + m2
 */
static ComplexMatrix add(ComplexMatrix m1, ComplexMatrix m2) {
    ComplexMatrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            m.a[i][j] = c.add(m1.a[i][j], m2.a[i][j]);
    return m;
}

/**
 * Gives differences of two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 - m2
 */
static ComplexMatrix sub(ComplexMatrix m1, ComplexMatrix m2) {
    ComplexMatrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            m.a[i][j] = c.sub(m1.a[i][j], m2.a[i][j]);
    return m;
}

/**
 * Gives product of two matrices
 * @param m1 first matrix
 * @param m2 second matrix
 * @return m1 x m2
 */
static ComplexMatrix dot(ComplexMatrix m1, ComplexMatrix m2) {
    ComplexMatrix m = null(m1.row, m2.col);
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            for (int k = 0; k < m1.col; k++)
                m.a[i][j] = c.add(m.a[i][j], c.dot(m1.a[i][k], m2.a[k][j]));
    return m;
}

/**
 * Gives scaled matrix
 * @param m matrix
 * @param s scale factor
 * @return s x m
 */
static ComplexMatrix scale(ComplexMatrix m, Complex s) {
    ComplexMatrix matrix = m;
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            matrix.a[i][j] = c.dot(matrix.a[i][j], s);
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
static ComplexMatrix minor(ComplexMatrix m, int i, int j) {
    ComplexMatrix matrix = null(m.row - 1, m.col - 1);
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
static Complex det(ComplexMatrix m) {
    if (m.row == 1)
        return m.a[0][0];
    Complex sum = {0, 0};
    for (int j = 0; j < m.col; j++)
        sum = c.add(sum,
                    c.dot(
                            c.dot(m.a[0][j], det(minor(m, 0, j))),
                            c.pow(c.get(-1.0f, 0.0f), c.get((float) j, 0.0f))
                    )
        );
    return sum;
}

/**
 * Gives cofactor of matrix
 * @param m matrix
 * @param i row index
 * @param j column index
 * @return cofactor of m at (i,j)
 */
static Complex cofactor(ComplexMatrix m, int i, int j) {
    return c.dot(
            c.pow(c.get(-1.0f, 0.0f), c.get((float) (i + j), 0.0f)),
            det(minor(m, i, j))
    );
}

/**
 * Gives cofactor matrix
 * @param m given matrix
 * @return cofactor matrix of m
 */
static ComplexMatrix cofactorMatrix(ComplexMatrix m) {
    ComplexMatrix matrix = null(m.row, m.col);
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
static ComplexMatrix transpose(ComplexMatrix m) {
    ComplexMatrix matrix = null(m.col, m.row);
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
static ComplexMatrix adjoint(ComplexMatrix m) {
    return transpose(cofactorMatrix(m));
}

/**
 * Gives inverse of matrix
 * @param m given matrix
 * @return adjoint(m)/det(m)
 */
static ComplexMatrix inverse(ComplexMatrix m) {
    if (m.row == 1 && m.col == 1) {
        ComplexMatrix matrix = {.row=1, .col=1, c.div(c.get(1, 0), m.a[0][0])};
        return matrix;
    }
    return scale(adjoint(m), c.div(c.get(1, 0), det(m)));
}

/**
 * Gives normalize matrix
 * @param m given matrix
 * @return normalize each column or vector in matrix m
 */
static ComplexMatrix normalise(ComplexMatrix m) {
    ComplexMatrix matrix = m;
    for (int j = 0; j < matrix.col; j++) {
        Complex sum = c.get(0, 0);
        for (int i = 0; i < matrix.row; i++)
            sum = c.add(sum, c.dot(matrix.a[i][j], c.conjugate(matrix.a[i][j])));
        sum = c.pow(sum, c.getFromReal(0.5f));
        for (int i = 0; i < matrix.row; i++)
            matrix.a[i][j] = c.div(matrix.a[i][j], sum);
    }
    return matrix;
}

/**
 * Gives conjugate matrix
 * @param m given matrix
 * @return conjugate each element and transpose the resultant matrix
 */
static ComplexMatrix conjugate(ComplexMatrix m) {
    ComplexMatrix matrix = {.row=m.col,.col=m.row};
    for (int i = 0; i < matrix.row; i++)
        for (int j = 0; j < matrix.col; j++)
            matrix.a[i][j] = c.conjugate(m.a[j][i]);
    return matrix;
}

/**
 * Gives rank of matrix
 * @param m given matrix
 * @return rank of m
 */
static int rank(ComplexMatrix m) {
    ComplexMatrix A;
    if (m.row <= m.col)
        A = dot(m, transpose(m));
    else
        A = dot(transpose(m), m);
    A = normalise(A);
    for (int rank = A.row; rank > 0; rank--) {
        if (c.mag(det(A)) < 0.001f) {
            ComplexMatrix temp = identity(A.row - 1);
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
static ComplexMatrix givenRotationMatrix(ComplexMatrix m,int i1, int i2) {
    int j = i1;
    ComplexMatrix matrix = identity(m.row);
    Complex r = c.pow(
            c.add(
                    c.dot(m.a[i1][j], m.a[i1][j]),
                    c.dot(m.a[i2][j], m.a[i2][j])
            ),
            c.get(0.5f, 0)
    );
    int diff = i2 - i1;
    Complex cos = c.div(m.a[i1][j], r);
    Complex sine = c.div(m.a[i2][j], r);
    matrix.a[i1][j] = cos;
    matrix.a[i1][j + diff] = sine;
    matrix.a[i2][j] = c.dot(c.getFromReal(-1), sine);
    matrix.a[i2][j + diff] = cos;
    return matrix;
}

/**
 * This QR decomposition using Given rotation matrix (m = Q x R)
 * @param m given matrix
 * @param Q Orthogonal matrix
 * @param R Right triangular matrix
 */
static void getQRWithGiven(ComplexMatrix m, ComplexMatrix *Q, ComplexMatrix *R) {
    *R = m;
    *Q = identity(m.row);
    for (int j = 0; j < m.col; j++) {
        for (int i = 0; i < m.col; i++) {
            if (i <= j)
                continue;
            ComplexMatrix G = givenRotationMatrix(*R, j, i);
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
static void getQRWithGram(ComplexMatrix m, ComplexMatrix *Q, ComplexMatrix *R) {
    //TODO not working with complex number
    *Q = null(m.row,m.row);
    *R = null(m.row,m.row);

    ComplexMatrix mt = conjugate(m);
    ComplexMatrix a_t[mt.col];
    ComplexMatrix u_n[mt.col];
    for (int i = 0; i < mt.row; i++) {
        a_t[i] = null(1,mt.col);
        for (int j = 0; j < mt.col; ++j)
            a_t[i].a[0][j]=mt.a[i][j];
        u_n[i] = null(mt.row, 1);
    }

    for (int k = 0; k < mt.col; k++) {
        ComplexMatrix proj = scale(u_n[0],dot(a_t[k],u_n[0]).a[0][0]);
        u_n[k]=sub(conjugate(a_t[k]),proj);
        for (int t = 1; t < k-1; ++t) {
            proj = scale(u_n[t],dot(u_n[k],u_n[t]).a[0][0]);
            u_n[k] = sub(conjugate(u_n[k]),proj);
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
static int isIdentityMatrix(ComplexMatrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if ((i == j) &&
                (c.mag(c.sub(m.a[i][j], c.get(1, 0))) > COMPLEX_EPSILON)) {
                return 0;
            } else if (c.mag(m.a[i][j]) > COMPLEX_EPSILON)
                return 0;
    return 1;
}

/**
 * Check if given matrix is diagonal matrix
 * @param m given matrix
 * @return returns 1 if m is diagonal matrix
 */
static int isDiagonalMatrix(ComplexMatrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if ((i != j) && (c.mag(m.a[i][j]) > COMPLEX_EPSILON))
                return 0;
    return 1;
}

/**
 * Check if given matrix is null matrix
 * @param m given matrix
 * @return returns 1 if m is null matrix
 */
static int isNullMatrix(ComplexMatrix m) {
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            if (c.mag(m.a[i][j]) > COMPLEX_EPSILON)
                return 0;
    return 1;
}


static Complex eigenCharacteristicsFunction(Complex x, void **args) {
    ComplexMatrix *m = (ComplexMatrix *) args[0];
    return det(add(*m, scale(identity(m->row), x)));
}

/**
 * Gives eigen decomposition of given matrix (m = V x D)
 * @param m given matrix
 * @param V eigen vector matrix
 * @param D eigen value matrix
 */
static void getEigenDecompose(ComplexMatrix m, ComplexMatrix *V, ComplexMatrix *D) {
    if (isDiagonalMatrix(m)) {
        *D = m;
        *V = identity(m.row);
        return;
    }

    //Find magnitude of each vector
    Complex mag = c.get(0, 0);
    for (int j = 0; j < m.col; j++)
        for (int i = 0; i < m.row; i++)
            mag = c.add(mag, c.dot(m.a[i][j], c.conjugate(m.a[i][j])));
    mag = c.pow(mag, c.getFromReal(0.5f));

    // Scale matrix
    m = scale(m, c.div(c.getFromReal(1.0f), mag));

    Complex roots[m.row];
    void *args[1] = {&m};
    StaticPolynomial.getRoots(m.row, eigenCharacteristicsFunction, args, roots);
    for (int i = 0; i < m.row; ++i) {
        for (int j = i + 1; j < m.row; ++j) {
            if (c.mag(roots[i]) > c.mag(roots[j]))
                continue;
            Complex temp = roots[j];
            roots[j] = roots[i];
            roots[i] = temp;
        }
    }

    *D = null(m.row, m.col);
    *V = identity(m.row);
    for (int i = 0; i < D->row; i++)
        D->a[i][i] = c.sub(c.getFromReal(0), roots[i]);

    // Re-scaling
    *D = scale(*D, mag);
    m = scale(m, mag);

    for (int n = 0; n < D->row; n++) {
        ComplexMatrix A = sub(m, scale(identity(D->row), D->a[n][n]));
        ComplexMatrix X;
        if (!isNullMatrix(A)) {
            X = StaticLinear.solveHomogeneous(A);
        } else {
            X = null(A.row, 1);
            X.a[n][0] = c.get(1.0f, 0.0f);
        }
        X = normalise(X);
        for (int i = 0; i < V->col; ++i)
            V->a[i][n] = X.a[i][0];
    }
}

/**
 * Gives singular value decomposition of given matrix (m = U * S x conjugate(V))
 * @param m given matrix of size (r x c)
 * @param U Unitary matrix of size (r x r)
 * @param S Rectangular diagonal matrix of size (r x c)
 * @param V Unitary matrix of size (c x c)
 */
static void getSVDDecompose(ComplexMatrix m, ComplexMatrix *U, ComplexMatrix *S, ComplexMatrix *V) {
    ComplexMatrix A = m;
    if (A.row >= A.col) {
        ComplexMatrix AAt = dot(A, conjugate(A));
        StaticComplexMatrix.print(AAt);
        printf("\n");

        ////////////LEFT SINGULAR MATRIX//////////////////
        ComplexMatrix D = AAt;
        ComplexMatrix Q, R;
        *U = identity(D.row);
        for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; ++i) {
            getQRWithGiven(D, &Q, &R);
            D = dot(R, Q);
            *U = dot(*U, Q);
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = c.pow(D.a[i][i], c.getFromReal(0.5f));

        ////////////RIGHT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < min(S->row, S->col); i++)
            S_inv.a[i][i] = c.div(c.getFromReal(1), S_inv.a[i][i]);
        *V = dot(S_inv, dot(conjugate(*U), A));
        *V = conjugate(*V);
    } else {
        ComplexMatrix AtA = dot(conjugate(A), A);
        ////////////RIGHT SINGULAR MATRIX//////////////////
        ComplexMatrix D = AtA;
        ComplexMatrix Q, R;
        *V = identity(D.row);
        for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; ++i) {
            getQRWithGiven(D, &Q, &R);
            D = dot(R, Q);
            *V = dot(*V, Q);
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = c.pow(D.a[i][i], c.get(0.5f, 0));

        ////////////LEFT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < S_inv.col; i++)
            S_inv.a[i][i] = c.div(c.get(1, 0), S_inv.a[i][i]);
        *U = dot(A, dot(*V, S_inv));
    }
}

///////////////////////////OPERATIONS///////////////////////////////////////////////////
/**
 * Returns matrix applying function f to each diagonal element
 * @param f function that take complex number and return complex number
 * @param m give matrix
 * @return modified matrix
 */
static ComplexMatrix diagonalOperation(Complex (*f)(Complex x), ComplexMatrix m) {
    for (int i = 0; i < min(m.row, m.col); ++i)
        m.a[i][i] = f(m.a[i][i]);
    return m;
}

/**
 * Returns matrix applying function f to each element
 * @param f function that takes complex number and return complex number
 * @param m give matrix
 * @return modified matrix
 */
static ComplexMatrix elementWiseOperation(Complex (*f)(Complex x), ComplexMatrix m) {
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            m.a[i][j] = f(m.a[i][j]);
    return m;
}

/**
 * This returns matrix after applying analytical function (should be taylor expansion!!!)
 * @param f analytic function that takes complex number and return complex number
 * @param m give matrix
 * @return modified matrix
 */
static ComplexMatrix analyticFunction(Complex (*f)(Complex x), ComplexMatrix m) {
    if (isDiagonalMatrix(m))
        return diagonalOperation(f, m);

    ComplexMatrix V, D;
    getEigenDecompose(m, &V, &D);
    return dot(dot(V, diagonalOperation(f, D)), inverse(V));
}

/**
 * This returns matrix after applying analytical function (should be taylor expansion!!!)
 * @param f analytic function that takes two complex number and return complex number
 * @param m given matrix
 * @param c1 given complex number
 * @return modified matrix
 */
static ComplexMatrix analyticFunction2(Complex (*f)(Complex x, Complex c1), ComplexMatrix m, Complex c1) {
    if (isDiagonalMatrix(m)) {
        for (int i = 0; i < min(m.row, m.col); ++i)
            m.a[i][i] = f(m.a[i][i], c1);
        return m;
    }

    ComplexMatrix V, D;
    getEigenDecompose(m, &V, &D);
    for (int i = 0; i < min(D.row, D.col); ++i)
        D.a[i][i] = f(D.a[i][i], c1);
    return dot(dot(V, D), inverse(V));
}

//////////////////////////EXPONENTIAL, LOGARITHMIC AND POWER/////////////////////////////
/**
 * Gives exponential of matrix
 * @param m given matrix
 * @return exp(m)
 */
static ComplexMatrix expm(ComplexMatrix m) {
    return analyticFunction(c.exp, m);
}

/**
 * Gives natural log of matrix
 * @param m given matrix
 * @return ln(m)
 */
static ComplexMatrix lnm(ComplexMatrix m) {
    return analyticFunction(c.ln, m);
}

/**
 * Gives pow of matrix raise to complex number
 * @param m given matrix
 * @param c given complex number
 * @return m^c
 */
static ComplexMatrix powm(ComplexMatrix m, Complex c1) {
    return analyticFunction2(c.pow, m, c1);
}

////////////////////////PRINTING////////////////////////////////////
/**
 * Print the matrix in form of x + iy
 * @param m given matrix
 */
static void print(ComplexMatrix m) {
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j) {
            c.print(m.a[i][j]);
            printf("    ");
        }
        printf("\n");
    }
}

/**
 * Print the matrix in form of r < theta
 * @param m given matrix
 */
static void printA(ComplexMatrix m) {
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j) {
            c.printA(m.a[i][j]);
            printf("    ");
        }
        printf("\n");
    }
}

/////////////////////////CONVERSION///////////////////
/**
 * Convert complex matrix to real matrix
 * @param m given complex matrix
 * @return only magnitude of complex matrix
 */
static Matrix toMatrix(ComplexMatrix m) {
    Matrix matrix = {.row = m.row, .col=m.col};
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            matrix.a[i][j] = m.a[i][j].x;
    return matrix;
}

/**
 * Convert complex matrix to real matrix
 * @param m given complex matrix
 * @return real part of complex matrix
 */
static Matrix toMatrixX(ComplexMatrix m) {
    Matrix matrix = {.row = m.row, .col=m.col};
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            matrix.a[i][j] = c.mag(m.a[i][j]);
    return matrix;
}

/**
 * Convert complex matrix to real matrix
 * @param m given complex matrix
 * @return imaginary part of complex matrix
 */
static Matrix toMatrixY(ComplexMatrix m) {
    Matrix matrix = {.row = m.row, .col=m.col};
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            matrix.a[i][j] = m.a[i][j].y;
    return matrix;
}

/**
 * Convert complex matrix to real matrix
 * @param m given complex matrix
 * @return magnitude part of complex matrix
 */
static Matrix toMatrixMag(ComplexMatrix m) {
    Matrix matrix = {.row = m.row, .col=m.col};
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            matrix.a[i][j] = c.mag(m.a[i][j]);
    return matrix;
}

/**
 * Convert complex matrix to real matrix
 * @param m given complex matrix
 * @return phase part of complex matrix
 */
static Matrix toMatrixPhase(ComplexMatrix m) {
    Matrix matrix = {.row = m.row, .col=m.col};
    for (int i = 0; i < m.row; ++i)
        for (int j = 0; j < m.col; ++j)
            matrix.a[i][j] = c.angle(m.a[i][j]);
    return matrix;
}

__attribute__((unused)) struct ComplexMatrixControl StaticComplexMatrix = {
        .get = get,
        .getFromMatrix = getFromMatrix,
        .getFromReal = getFromReal,
        .getXY = getXY,
        .getConstantMatrix = getConstantMatrix,
        .getConstantMatrixFromReal = getConstantMatrixFromReal,
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
        .conjugate = conjugate,
        .rank = rank,
        .givenRotationMatrix = givenRotationMatrix,
        .getQRWithGiven = getQRWithGiven,
        .getQRWithGram = getQRWithGram,
        .isIdentityMatrix = isIdentityMatrix,
        .isDiagonalMatrix = isDiagonalMatrix,
        .isNullMatrix = isNullMatrix,
        .getEigenDecompose = getEigenDecompose,
        .getSVDDecompose = getSVDDecompose,
        .diagonalOperation = diagonalOperation,
        .elementWiseOperation = elementWiseOperation,
        .analyticFunction = analyticFunction,
        .analyticFunction2 = analyticFunction2,
        .exp = expm,
        .ln = lnm,
        .pow = powm,
        .print = print,
        .printA = printA,
        .toMatrix = toMatrix,
        .toMatrixX = toMatrixX,
        .toMatrixY = toMatrixY,
        .toMatrixMag = toMatrixMag,
        .toMatrixPhase = toMatrixPhase
};
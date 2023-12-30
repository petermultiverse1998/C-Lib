//
// Created by peter on 12/29/2023.
//

#include <stdio.h>
#include "complex_matrix.h"
#include "polynomial.h"
#include "linear.h"

#define min(a, b) (a<=b?a:b)

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
 * @param row row of matrix
 * @param col column of matrix
 * @param a element of matrix in 2D array of real number
 * @return complex matrix
 */
static ComplexMatrix getFromReal(int row, int col, float a[row][col]) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = StaticComplex.getFromReal(a[i][j]);
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
            matrix.a[i][j] = StaticComplex.get(x[i][j], y[i][j]);
    return matrix;
}

/**
 * This will give complex matrix
 * @param row row of matrix
 * @param col column of matrix
 * @param c complex number for all the elements in matrix
 * @return complex matrix x+iy
 */
static ComplexMatrix getConstantMatrix(int row, int col, Complex c) {
    ComplexMatrix matrix = {.row=row, .col=col};
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            matrix.a[i][j] = c;
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
            matrix.a[i][j] = StaticComplex.getFromReal(x);
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
            m.a[i][j] = StaticComplex.add(m1.a[i][j], m2.a[i][j]);
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
            m.a[i][j] = StaticComplex.sub(m1.a[i][j], m2.a[i][j]);
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
                m.a[i][j] = StaticComplex.add(m.a[i][j], StaticComplex.dot(m1.a[i][k], m2.a[k][j]));
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
            matrix.a[i][j] = StaticComplex.dot(matrix.a[i][j], s);
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
        sum = StaticComplex.add(sum,
                                StaticComplex.dot(
                                        StaticComplex.dot(m.a[0][j], det(minor(m, 0, j))),
                                        StaticComplex.pow(StaticComplex.get(-1.0f, 0.0f),
                                                          StaticComplex.get((float) j, 0.0f))
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
    return StaticComplex.dot(
            StaticComplex.pow(StaticComplex.get(-1.0f, 0.0f), StaticComplex.get((float) (i + j), 0.0f)),
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
        ComplexMatrix matrix = {.row=1, .col=1, StaticComplex.div(StaticComplex.get(1, 0), m.a[0][0])};
        return matrix;
    }
    return scale(adjoint(m), StaticComplex.div(StaticComplex.get(1, 0), det(m)));
}

/**
 * Gives normalize matrix
 * @param m given matrix
 * @return normalize each column or vector in matrix m
 */
static ComplexMatrix normalise(ComplexMatrix m) {
    ComplexMatrix matrix = m;
    for (int j = 0; j < matrix.col; j++) {
        Complex sum = StaticComplex.get(0, 0);
        for (int i = 0; i < matrix.row; i++)
            sum = StaticComplex.add(sum, StaticComplex.dot(matrix.a[i][j], StaticComplex.conjugate(matrix.a[i][j])));
        sum = StaticComplex.pow(sum,StaticComplex.getFromReal(0.5f));
        for (int i = 0; i < matrix.row; i++)
            matrix.a[i][j] = StaticComplex.div(matrix.a[i][j], sum);
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
            matrix.a[i][j] = StaticComplex.conjugate(m.a[j][i]);
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
        if (StaticComplex.mag(det(A)) < 0.001f) {
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
    Complex r = StaticComplex.pow(
            StaticComplex.add(
                    StaticComplex.pow(m.a[i1][j], StaticComplex.get(2, 0)),
                    StaticComplex.pow(m.a[i2][j], StaticComplex.get(2, 0))
            ),
            StaticComplex.get(0.5f, 0)
    );
    int diff = i2 - i1;
    Complex c = StaticComplex.div(m.a[i1][j], r);
    Complex s = StaticComplex.div(StaticComplex.dot(StaticComplex.get(-1, 0), m.a[i2][j]), r);
    matrix.a[i1][j] = c;
    matrix.a[i1][j + diff] = StaticComplex.dot(StaticComplex.get(-1, 0), s);
    matrix.a[i2][j] = s;
    matrix.a[i2][j + diff] = c;
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
                (StaticComplex.mag(StaticComplex.sub(m.a[i][j], StaticComplex.get(1, 0))) > COMPLEX_EPSILON)) {
                return 0;
            } else if (StaticComplex.mag(m.a[i][j]) > COMPLEX_EPSILON)
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
            if ((i != j) && (StaticComplex.mag(m.a[i][j]) > COMPLEX_EPSILON))
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
            if (StaticComplex.mag(m.a[i][j]) > COMPLEX_EPSILON)
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
    Complex roots[m.row];
    void *args[1] = {&m};
    StaticPolynomial.getRoots(m.row, eigenCharacteristicsFunction, args, roots);

    *D = null(m.row, m.col);
    *V = identity(m.row);
    for (int i = 0; i < D->row; i++)
        D->a[i][i] = StaticComplex.sub(StaticComplex.get(0, 0), roots[i]);
    for (int n = 0; n < D->row; n++) {
        ComplexMatrix A = sub(m, scale(identity(D->row), D->a[n][n]));
        ComplexMatrix X;
        if (!isNullMatrix(A)) {
            X = StaticLinear.solveHomogeneous(A);
        } else {
            X = null(A.row, 1);
            X.a[n][0] = StaticComplex.get(1.0f, 0.0f);
        }
        X = normalise(X);
        for (int i = 0; i < V->col; ++i)
            V->a[i][n] = X.a[i][0];
    }
}

/**
 * Gives eigen decomposition of given matrix using QR method(m = V x D).
 * Only works for symmetric matrix
 * @param m given matrix
 * @param V eigen vector matrix
 * @param D eigen value matrix
 */
static void getEigenDecomposeQR(ComplexMatrix m, ComplexMatrix *V, ComplexMatrix *D) {
    if (isDiagonalMatrix(m)) {
        *D = m;
        *V = identity(m.row);
        return;
    }
    *D = m;
    ComplexMatrix Q;
    ComplexMatrix R;
    *V = identity(m.row);
    for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; i++) {
        getQRWithGiven(*D, &Q, &R);
        *D = dot(R, Q);
    }
    for (int i = 0; i < D->row; ++i)
        for (int j = 0; j < D->col; ++j)
            if (i != j)
                D->a[i][j] = StaticComplex.get(0, 0);

    for (int n = 0; n < D->row; n++) {
        ComplexMatrix A = sub(m, scale(identity(D->row), D->a[n][n]));
        ComplexMatrix X;
        if (!isNullMatrix(A)) {
            X = StaticLinear.solveHomogeneous(A);
        } else {
            X = null(A.row, 1);
            X.a[n][0] = StaticComplex.get(1.0f, 0.0f);
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

        ComplexMatrix D;
        getEigenDecompose(AAt, U, &D);
        ////////////LEFT SINGULAR MATRIX//////////////////

        //ASCENDING
        for (int i = 0; i < D.row - 1; i++) {
            for (int j = i; j < D.row - 1; j++) {
                if (StaticComplex.mag(D.a[j][j]) < StaticComplex.mag(D.a[j + 1][j + 1])) {
                    Complex temp = D.a[j][j];
                    D.a[j][j] = D.a[j + 1][j + 1];
                    D.a[j + 1][j + 1] = temp;

                    for (int k = 0; k < U->col; ++k) {
                        Complex tempV = U->a[j][k];
                        U->a[j][k] = U->a[j + 1][k];
                        U->a[j + 1][k] = tempV;
                    }
                }
            }
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = StaticComplex.pow(D.a[i][i], StaticComplex.get(0.5f, 0));

        ////////////RIGHT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < min(S->row, S->col); i++)
            S_inv.a[i][i] = StaticComplex.div(StaticComplex.get(1, 0), S->a[i][i]);
        *V = dot(S_inv, dot(conjugate(*U), A));
        *V = conjugate(*V);
    } else {
        ComplexMatrix AtA = dot(conjugate(A), A);
        ComplexMatrix D;
        getEigenDecompose(AtA, V, &D);
        ////////////RIGHT SINGULAR MATRIX//////////////////

        //ASCENDING
        for (int i = 0; i < D.row - 1; i++) {
            for (int j = i; j < D.row - 1; j++) {
                if (StaticComplex.mag(D.a[j][j]) < StaticComplex.mag(D.a[j + 1][j + 1])) {
                    Complex temp = D.a[j][j];
                    D.a[j][j] = D.a[j + 1][j + 1];
                    D.a[j + 1][j + 1] = temp;

                    for (int k = 0; k < U->col; ++k) {
                        Complex tempV = U->a[j][k];
                        V->a[j][k] = V->a[j + 1][k];
                        V->a[j + 1][k] = tempV;
                    }
                }
            }
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = StaticComplex.pow(D.a[i][i], StaticComplex.get(0.5f, 0));

        ////////////LEFT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < S_inv.col; i++)
            S_inv.a[i][i] = StaticComplex.div(StaticComplex.get(1, 0), S_inv.a[i][i]);
        *U = dot(A, dot(*V, S_inv));
    }
}

/**
 * Gives singular value decomposition of given matrix using QR method (m = U * S x conjugate(V))
 * @param m given matrix of size (r x c)
 * @param U Unitary matrix of size (r x r)
 * @param S Rectangular diagonal matrix of size (r x c)
 * @param V Unitary matrix of size (c x c)
 */
static void getSVDDecomposeQR(ComplexMatrix m, ComplexMatrix *U, ComplexMatrix *S, ComplexMatrix *V) {
    ComplexMatrix A = m;
    if (A.row >= A.col) {
        ComplexMatrix AAt = dot(A, conjugate(A));
        ComplexMatrix D = AAt;
        for (int i = 0; i < COMPLEX_MATRIX_QR_EIGEN_ITERATION; i++) {
            ComplexMatrix Q, R;
            getQRWithGiven(D, &Q, &R);
            D = dot(R, Q);
        }

        //ASCENDING
        for (int i = 0; i < D.row - 1; i++) {
            for (int j = i; j < D.row - 1; j++) {
                if (StaticComplex.mag(D.a[j][j]) < StaticComplex.mag(D.a[j + 1][j + 1])) {
                    Complex temp = D.a[j][j];
                    D.a[j][j] = D.a[j + 1][j + 1];
                    D.a[j + 1][j + 1] = temp;
                }
            }
        }

        ////////////LEFT SINGULAR MATRIX//////////////////
        *U = identity(AAt.row);
        ComplexMatrix temp;
        for (int n = 0; n < AAt.row; n++) {
            temp = sub(AAt, scale(identity(AAt.row), D.a[n][n]));
            ComplexMatrix X;
            if (!isNullMatrix(temp)) {
                X = StaticLinear.solveHomogeneous(temp);
            } else {
                X = null(A.row, 1);
                X.a[n][0] = StaticComplex.get(1.0f, 0.0f);
            }
            X = normalise(X);
            for (int j = 0; j < U->col; ++j)
                U->a[n][j] = X.a[j][0];
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = StaticComplex.pow(D.a[i][i], StaticComplex.get(0.5f, 0));

        ////////////RIGHT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < S_inv.row; i++)
            S_inv.a[i][i] = StaticComplex.div(StaticComplex.get(1, 0), S_inv.a[i][i]);
        *V = dot(S_inv, dot(conjugate(*U), A));
        *V = conjugate(*V);

    } else {
        ComplexMatrix AtA = dot(conjugate(A), A);
        ComplexMatrix D = AtA;
        for (int i = 0; i < 100; i++) {
            ComplexMatrix Q, R;
            getQRWithGiven(D, &Q, &R);
            D = dot(R, Q);
        }

        //ASCENDING
        for (int i = 0; i < D.row - 1; i++) {
            for (int j = i; j < D.row - 1; j++) {
                if (StaticComplex.mag(D.a[j][j]) < StaticComplex.mag(D.a[j + 1][j + 1])) {
                    Complex temp = D.a[j][j];
                    D.a[j][j] = D.a[j + 1][j + 1];
                    D.a[j + 1][j + 1] = temp;
                }
            }
        }

        ////////////RIGHT SINGULAR MATRIX//////////////////
        *V = identity(AtA.row);
        ComplexMatrix temp;
        for (int n = 0; n < AtA.row; n++) {
            temp = sub(AtA, scale(identity(AtA.row), D.a[n][n]));
            ComplexMatrix X;
            if (!isNullMatrix(temp)) {
                X = StaticLinear.solveHomogeneous(temp);
            } else {
                X = null(A.row, 1);
                X.a[n][0] = StaticComplex.get(1.0f, 0.0f);
            }
            for (int j = 0; j < V->col; ++j)
                V->a[n][j] = X.a[j][0];
        }

        ////////SINGULAR MATRIX//////////////
        *S = null(A.row, A.col);
        for (int i = 0; i < min(S->row, S->col); i++)
            S->a[i][i] = StaticComplex.pow(D.a[i][i], StaticComplex.get(0.5f, 0));

        ////////////RIGHT SINGULAR MATRIX//////////////////
        ComplexMatrix S_inv = transpose(*S);
        for (int i = 0; i < S_inv.col; i++)
            S_inv.a[i][i] = StaticComplex.div(StaticComplex.get(1, 0), S_inv.a[i][i]);
        *U = dot(S_inv, dot(*V, A));
    }
}

////////////////////////PRINTING////////////////////////////////////
/**
 * Print the matrix in form of x + iy
 * @param m given matrix
 */
static void print(ComplexMatrix m) {
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j) {
            StaticComplex.print(m.a[i][j]);
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
            StaticComplex.printA(m.a[i][j]);
            printf("    ");
        }
        printf("\n");
    }
}

__attribute__((unused)) struct ComplexMatrixControl StaticComplexMatrix = {
        .get = get,
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
        .getEigenDecomposeQR = getEigenDecomposeQR,
        .getSVDDecompose = getSVDDecompose,
        .getSVDDecomposeQR = getSVDDecomposeQR,
        .print = print,
        .printA = printA
};
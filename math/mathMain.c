//
// Created by peter on 12/27/2023.
//
#include <math.h>
#include "stdio.h"
#include "complex.h"
#include "polynomial.h"
#include "complex_matrix.h"
#include "linear.h"

#define END ;printf("\n")
#define c StaticComplex

void complexTest() {
    //Basic operations test
    printf("Basic Operations\n");
    c.print(c.get(1, 2))END;
    c.print(c.getFromReal(3))END;
    c.printA(c.getFromPolar(4, 1))END;

    //Complex operations test
    printf("\nComplex Operations\n");
    Complex c1 = c.get(1, 3);
    c.print(c1)END;
    c.print(c.conjugate(c1))END;
    printf("%f\n",c.mag(c1));
    printf("%f\n",c.angle(c1));
    c.print(c.normalise(c1))END;

    //Arithmetic operations test
    printf("\nArithmetic Operations\n");
    c1 = c.get(2,5);
    Complex c2 = c.get(3,4);
    c.print(c1)END;
    c.print(c2)END;
    c.print(c.add(c1,c2))END;
    c.print(c.sub(c1,c2))END;
    c.print(c.dot(c1,c2))END;
    c.print(c.div(c1,c2))END;

    //Exponential and power operations test
    printf("\nExponential and Power Operations\n");
    c1 = c.get(2,7);
    c2 = c.get(4,1);
    c.print(c1)END;
    c.print(c2)END;
    c.print(c.exp(c1))END;
    c.print(c.ln(c1))END;
    c.print(c.pow(c1,c2))END;

    //Trigonometric operations test
    printf("\nTrigonometric Operations\n");
    c1 = c.get(2,7);
    c.print(c1)END;
    c.print(c.sin(c1))END;
    c.print(c.cos(c1))END;
    c.print(c.tan(c1))END;

    //Inverse Trigonometric operations test
    printf("\nInverse Trigonometric Operations\n");
    c1 = c.get(2,7);
    c2 = c.get(4,1);
    c.print(c1)END;
    c.print(c2)END;
    c.print(c.asin(c1))END;
    c.print(c.acos(c1))END;
    c.print(c.atan(c1))END;
    c.print(c.atan2(c2,c1))END;

    //Hyperbolic Trigonometric operations test
    printf("\nHyperbolic Trigonometric Operations\n");
    c1 = c.get(2,7);
    c.print(c1)END;
    c.print(c.sinh(c1))END;
    c.print(c.cosh(c1))END;
    c.print(c.tanh(c1))END;

    //Inverse Hyperbolic Trigonometric operations test
    printf("\nInverse Hyperbolic Trigonometric Operations\n");
    c1 = c.get(2,7);
    c.print(c1)END;
    c.print(c.asinh(c1))END;
    c.print(c.acosh(c1))END;
    c.print(c.atanh(c1))END;
}

#define degree 8
Complex f(Complex x,void**args){
    Complex c1[degree+1];
    for (int i = 0; i < (degree+1); ++i)
        c1[i] = c.getFromReal(0);
    c1[0] = c.getFromReal(-1);
    c1[degree] = c.getFromReal(1);


    Complex sum = c.get(0,0);
    for (int i = 0; i < (degree+1); ++i)
        sum = c.add(sum,c.dot(c1[i],c.pow(x,c.getFromReal((float)i))));
    return sum;
}
void polynomialTest(){
    Complex r[degree];
    StaticPolynomial.getRoots(degree,f,NULL,r);
    for(int i=0;i<degree;i++) {
        c.printA(r[i])END;
    }
}

/**************COMPLEX MATRIX TEST**********/
Complex a[3][3] = {
        {{3, 5}, {1, 6}, {4, 2}},
        {{1, 9}, {5, 5}, {9, 1}},
        {{2, 4}, {6, 1}, {5, 3}}
};
float x[3][3] = {
        {5, 6, 2},
        {9, 5, 1},
        {4, 1, 3}
};
float y[3][3] = {
        {5, 6, 2},
        {9, 5, 1},
        {9, 5, 1}
};

#define cm StaticComplexMatrix
void complexMatrixTestBasic() {
    printf("\nBasic\n");
    cm.print(cm.get(3, 3, a))END;
    cm.print(cm.getFromReal(3, 3, x))END;
    cm.print(cm.getXY(3, 3, x, x))END;
    cm.print(cm.getConstantMatrix(1, 3, c.get(3, 1)))END;
    cm.print(cm.getConstantMatrixFromReal(3, 1, 3.14f))END;
}

void complexMatrixSpecialMatrixTest(){
    printf("\nSpecial Matrix\n");
    cm.print(cm.identity(3))END;
    cm.print(cm.null(3,3))END;
}

void complexMatrixArithmeticTest(){
    printf("\nArithmetic Operations\n");
    ComplexMatrix m1 = cm.get(3,3,a);
    ComplexMatrix m2 = cm.getFromReal(3,3,x);
    cm.print(m1)END;
    cm.print(m2)END;
    cm.print(cm.add(m1,m2))END;
    cm.print(cm.sub(m1,m2))END;
    cm.print(cm.dot(m1,m2))END;
    cm.print(cm.scale(m1,c.get(3,1)))END;
}

void complexMatrixOperationTest(){
    printf("\nMatrix Operations\n");
    ComplexMatrix m1 = cm.get(2,3,a);
//    ComplexMatrix m1 = cm.getFromReal(3,3,x);
//    ComplexMatrix m2 = cm.getFromReal(3,3,y);
    cm.print(m1)END;
    cm.print(cm.minor(m1,1,1))END;
    c.print(cm.det(m1))END;
    c.print(cm.cofactor(m1,1,1))END;
    cm.print(cm.cofactorMatrix(m1))END;
    cm.print(cm.transpose(m1))END;
    cm.print(cm.adjoint(m1))END;
    cm.print(cm.inverse(m1))END;
    cm.print(cm.normalise(m1))END;
    cm.print(cm.conjugate(m1))END;
    printf("%d\n",cm.rank(cm.getFromReal(3,3,y)));
}

void complexMatrixDecompositionTest(){
    printf("\nDecomposition\n");
//    ComplexMatrix m1 = cm.get(3,3,a);
    ComplexMatrix m1 = cm.getFromReal(3,3,x);
//    ComplexMatrix m2 = cm.getFromReal(3,3,y);

    printf("Given \n");
    ComplexMatrix G = cm.givenRotationMatrix(m1,1,2);
    cm.print(G)END;
    ComplexMatrix A = cm.dot(G,m1);
    cm.print(A)END;

    printf("QR with given\n");
    ComplexMatrix Q,R;
    cm.getQRWithGiven(m1,&Q,&R);
    cm.print(Q)END;
    cm.print(R)END;
    cm.print(cm.dot(Q,R))END;

    //TODO QR with Gram not working with complex number
    printf("QR with gram\n");
    cm.getQRWithGram(m1,&Q,&R);
    cm.print(Q)END;
    cm.print(R)END;
    cm.print(cm.dot(Q,R))END;

    printf("Eigen Decomposition\n");
    ComplexMatrix V,D;
    cm.getEigenDecompose(m1,&V,&D);
    cm.print(V)END;
    cm.print(D)END;
    cm.print(cm.dot(cm.dot(V,D),cm.inverse(V)))END;

    //TODO works only for real matrix
    printf("Singular Value Decomposition\n");
    ComplexMatrix U;
    cm.getSVDDecompose(m1,&U,&D,&V);
    cm.print(U)END;
    cm.print(D)END;
    cm.print(V)END;
    cm.print(cm.dot(cm.dot(U,D),cm.conjugate(V)))END;
}

void complexMatrixFunctionTest(){
    ComplexMatrix temp;

    printf("\nFunction Mapping\n");
    ComplexMatrix m1 = cm.get(3,3,a);
//    ComplexMatrix m1 = cm.getFromReal(3,3,x);
    cm.print(m1)END;
    m1 = cm.scale(m1,c.getFromReal(0.1f));
    cm.print(cm.diagonalOperation(c.exp,m1))END;
    cm.print(cm.elementWiseOperation(c.exp,m1))END;

    printf("\nAnalytic Function\n");
    ComplexMatrix V,D;
    temp = cm.analyticFunction(c.exp,m1);
    cm.print(temp)END;
    cm.print(cm.ln(temp))END;

    printf("\nAnalytic Function 2\n");
    temp = cm.analyticFunction2(c.pow,m1,c.getFromReal(0.5f));
    cm.print(temp)END;
    cm.print(cm.pow(temp,c.getFromReal(2.0f)))END;

    printf("\nExponential, Logarithmic and power\n");
    temp = cm.exp(m1);
    cm.print(temp)END;
    cm.print(cm.ln(temp))END;
    temp = cm.pow(m1,c.getFromReal(0.5f));
    cm.print(temp)END;
    cm.print(cm.pow(temp,c.getFromReal(2.0f)))END;
}

void complexMatrixTest() {
//    complexMatrixTestBasic();
//    complexMatrixSpecialMatrixTest();
//    complexMatrixArithmeticTest();
//    complexMatrixOperationTest();
//    complexMatrixDecompositionTest();
//    complexMatrixFunctionTest();
}

/**************LINEAR TEST**********/
void linearTest(){
    float x[2][2]={
            {5,9},
            {5,9}
    };
    ComplexMatrix A = cm.getFromReal(2,2,x);
    cm.print(A)END;
    ComplexMatrix X = StaticLinear.solveHomogeneous(A);
    cm.print(X);
}

/**************MATRIX TEST**********/
#define m StaticMatrix
void matrixTestBasic() {
    printf("\nBasic\n");
    m.print(m.get(3, 3, x))END;
}

void matrixSpecialMatrixTest(){
    printf("\nSpecial Matrix\n");
    m.print(m.identity(3))END;
    m.print(m.null(3,3))END;
}

void matrixArithmeticTest(){
    printf("\nArithmetic Operations\n");
    Matrix m1 = m.get(3,3,x);
    Matrix m2 = m.get(3,3,y);
    m.print(m1)END;
    m.print(m2)END;
    m.print(m.add(m1,m2))END;
    m.print(m.sub(m1,m2))END;
    m.print(m.dot(m1,m2))END;
    m.print(m.scale(m1,3))END;
}

void matrixOperationTest(){
    printf("\nMatrix Operations\n");
    Matrix m1 = m.get(3,3,x);
    m.print(m1)END;
    m.print(m.minor(m1,1,1))END;
    printf("%f\n",m.det(m1));
    printf("%f\n",m.cofactor(m1,1,1));
    m.print(m.cofactorMatrix(m1))END;
    m.print(m.transpose(m1))END;
    m.print(m.adjoint(m1))END;
    m.print(m.inverse(m1))END;
    m.print(m.normalise(m1))END;
    printf("%d\n",m.rank(m.get(3,3,y)));
}

void matrixDecompositionTest(){
    printf("\nDecomposition\n");
    Matrix m1 = m.get(3,3,x);

    printf("Given \n");
    Matrix G = m.givenRotationMatrix(m1,1,2);
    m.print(G)END;
    Matrix A = m.dot(G,m1);
    m.print(A)END;

    printf("QR with given\n");
    Matrix Q,R;
    m.getQRWithGiven(m1,&Q,&R);
    m.print(Q)END;
    m.print(R)END;
    m.print(m.dot(Q,R))END;

    //TODO QR with Gram not working with complex number
    printf("QR with gram\n");
    m.getQRWithGram(m1,&Q,&R);
    m.print(Q)END;
    m.print(R)END;
    m.print(m.dot(Q,R))END;

    //TODO works only for real matrix
    printf("Singular Value Decomposition\n");
    Matrix V,D,U;
    m.getSVDDecompose(m1,&U,&D,&V);
    m.print(U)END;
    m.print(D)END;
    m.print(V)END;
    m.print(m.dot(m.dot(U,D),m.transpose(V)))END;
}

void matrixFunctionTest(){
    Matrix temp;

    printf("\nFunction Mapping\n");
    Matrix m1 = m.get(3,3,x);
//    ComplexMatrix m1 = m.getFromReal(3,3,x);
    m1 = m.scale(m1,0.1f);
    m.print(m1)END;
    m.print(m.diagonalOperation(expf,m1))END;
    m.print(m.elementWiseOperation(expf,m1))END;

    printf("\nAnalytic Function\n");
    temp = m.analyticFunction(c.exp,m1);
    m.print(temp)END;
    m.print(m.ln(temp))END;

    printf("\nAnalytic Function 2\n");
    temp = m.analyticFunction2(c.pow,m1,0.5f);
    m.print(temp)END;
    m.print(m.pow(temp,2.0f))END;

    printf("\nExponential, Logarithmic and power\n");
    temp = m.exp(m1);
    m.print(temp)END;
    m.print(m.ln(temp))END;
    temp = m.pow(m1,0.5f);
    m.print(temp)END;
    m.print(m.pow(temp,2.0f))END;
}

void matrixTest() {
//    matrixTestBasic();
//    matrixSpecialMatrixTest();
//    matrixArithmeticTest();
//    matrixOperationTest();
//    matrixDecompositionTest();
//    matrixFunctionTest();
}

int main() {
//    complexTest();
//    polynomialTest();
//    complexMatrixTest();
    matrixTest();
//    linearTest();
}
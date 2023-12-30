//
// Created by peter on 12/27/2023.
//
#include "stdio.h"
#include "complex.h"
#include "polynomial.h"
#include "complex_matrix.h"
#include "linear.h"

#define END ;printf("\n")

void complexTest() {
    //Basic operations test
    printf("Basic Operations\n");
    StaticComplex.print(StaticComplex.get(1, 2))END;
    StaticComplex.print(StaticComplex.getFromReal(3))END;
    StaticComplex.printA(StaticComplex.getFromPolar(4, 1))END;

    //Complex operations test
    printf("\nComplex Operations\n");
    Complex c = StaticComplex.get(1, 3);
    StaticComplex.print(c)END;
    StaticComplex.print(StaticComplex.conjugate(c))END;
    printf("%f\n",StaticComplex.mag(c));
    printf("%f\n",StaticComplex.angle(c));
    StaticComplex.print(StaticComplex.normalise(c))END;

    //Arithmetic operations test
    printf("\nArithmetic Operations\n");
    Complex c1 = StaticComplex.get(2,5);
    Complex c2 = StaticComplex.get(3,4);
    StaticComplex.print(c1)END;
    StaticComplex.print(c2)END;
    StaticComplex.print(StaticComplex.add(c1,c2))END;
    StaticComplex.print(StaticComplex.sub(c1,c2))END;
    StaticComplex.print(StaticComplex.dot(c1,c2))END;
    StaticComplex.print(StaticComplex.div(c1,c2))END;

    //Exponential and power operations test
    printf("\nExponential and Power Operations\n");
    c1 = StaticComplex.get(2,7);
    c2 = StaticComplex.get(4,1);
    StaticComplex.print(c1)END;
    StaticComplex.print(c2)END;
    StaticComplex.print(StaticComplex.exp(c1))END;
    StaticComplex.print(StaticComplex.ln(c1))END;
    StaticComplex.print(StaticComplex.pow(c1,c2))END;

    //Trigonometric operations test
    printf("\nTrigonometric Operations\n");
    c = StaticComplex.get(2,7);
    StaticComplex.print(c)END;
    StaticComplex.print(StaticComplex.sin(c))END;
    StaticComplex.print(StaticComplex.cos(c))END;
    StaticComplex.print(StaticComplex.tan(c))END;

    //Inverse Trigonometric operations test
    printf("\nInverse Trigonometric Operations\n");
    c1 = StaticComplex.get(2,7);
    c2 = StaticComplex.get(4,1);
    StaticComplex.print(c1)END;
    StaticComplex.print(c2)END;
    StaticComplex.print(StaticComplex.asin(c1))END;
    StaticComplex.print(StaticComplex.acos(c1))END;
    StaticComplex.print(StaticComplex.atan(c1))END;
    StaticComplex.print(StaticComplex.atan2(c2,c1))END;

    //Hyperbolic Trigonometric operations test
    printf("\nHyperbolic Trigonometric Operations\n");
    c = StaticComplex.get(2,7);
    StaticComplex.print(c)END;
    StaticComplex.print(StaticComplex.sinh(c))END;
    StaticComplex.print(StaticComplex.cosh(c))END;
    StaticComplex.print(StaticComplex.tanh(c))END;

    //Inverse Hyperbolic Trigonometric operations test
    printf("\nInverse Hyperbolic Trigonometric Operations\n");
    c = StaticComplex.get(2,7);
    StaticComplex.print(c)END;
    StaticComplex.print(StaticComplex.asinh(c))END;
    StaticComplex.print(StaticComplex.acosh(c))END;
    StaticComplex.print(StaticComplex.atanh(c))END;
}

#define degree 8
Complex f(Complex x,void**args){
    Complex c[degree+1];
    for (int i = 0; i < (degree+1); ++i)
        c[i] = StaticComplex.getFromReal(0);
    c[0] = StaticComplex.getFromReal(-1);
    c[degree] = StaticComplex.getFromReal(1);


    Complex sum = StaticComplex.get(0,0);
    for (int i = 0; i < (degree+1); ++i)
        sum = StaticComplex.add(sum,StaticComplex.dot(c[i],StaticComplex.pow(x,StaticComplex.getFromReal((float)i))));
    return sum;
}
void polynomialTest(){
    Complex r[degree];
    StaticPolynomial.getRoots(degree,f,NULL,r);
    for(int i=0;i<degree;i++) {
        StaticComplex.printA(r[i])END;
    }
}

void matrixTest() {
    //Basic operations test
    Complex a[3][3]={
            {{3,5},{1,6},{4,2}},
            {{1,9},{5,5},{9,1}},
            {{2,4},{6,1},{5,3}}
    };

    float x[3][3]={
            {5,6,2},
            {9,5,1},
            {4,1,3}
    };
    float x1[3][2]={
            {5,6},
            {9,5},
            {4,1}
    };
    float y1[2][3]={
            {5,6,2},
            {9,5,1}
    };


    float y[3][3]={
            {5,6,2},
            {9,5,1},
            {9,5,1}
    };

//    printf("Basic\n");
//    StaticComplexMatrix.print(StaticComplexMatrix.get(3, 3, a))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.getFromReal(3, 3, x))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.getXY(3, 3, x,x))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.getConstantMatrix(1, 3, StaticComplex.get(3,1)))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.getConstantMatrixFromReal(3, 1, 3.14f))END;
//
//
//    printf("\nSpecial Matrix\n");
//    StaticComplexMatrix.print(StaticComplexMatrix.identity(3))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.null(3,3))END;

//    //Arithmetic operations test
//    printf("\nArithmetic Operations\n");
//    ComplexMatrix m1 = StaticComplexMatrix.get(3,3,a);
//    ComplexMatrix m2 = StaticComplexMatrix.getFromReal(3,3,x);
//    StaticComplexMatrix.print(m1)END;
//    StaticComplexMatrix.print(m2)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.add(m1,m2))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.sub(m1,m2))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(m1,m2))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.scale(m1,StaticComplex.get(3,1)))END;

    //Matrix Operations
    printf("\nMatrix Operations\n");
//    ComplexMatrix m = StaticComplexMatrix.get(2,3,a);
//    ComplexMatrix m = StaticComplexMatrix.getFromReal(3,2,x1);
    ComplexMatrix m = StaticComplexMatrix.getFromReal(2,3,y1);
//    ComplexMatrix m = StaticComplexMatrix.getFromReal(3,3,x);
    StaticComplexMatrix.print(m)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.minor(m,1,1))END;
//    StaticComplex.print(StaticComplexMatrix.det(m))END;
//    StaticComplex.print(StaticComplexMatrix.cofactor(m,1,1))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.cofactorMatrix(m))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.transpose(m))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.adjoint(m))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.inverse(m))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.normalise(m))END;
//    StaticComplexMatrix.print(StaticComplexMatrix.conjugate(m))END;
//    printf("%d\n",StaticComplexMatrix.rank(StaticComplexMatrix.getFromReal(3,3,y)));
//    ComplexMatrix G = StaticComplexMatrix.givenRotationMatrix(m,1,2);
//    StaticComplexMatrix.print(G)END;
//    ComplexMatrix A = StaticComplexMatrix.dot(G,m);
//    StaticComplexMatrix.print(A)END;
//    ComplexMatrix Q,R;
//    StaticComplexMatrix.getQRWithGiven(m,&Q,&R);
//    StaticComplexMatrix.print(Q)END;
//    StaticComplexMatrix.print(R)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(Q,R))END;
    //TODO QR with Gram not working with complex number
//    StaticComplexMatrix.getQRWithGram(m,&Q,&R);
//    StaticComplexMatrix.print(Q)END;
//    StaticComplexMatrix.print(R)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(Q,R))END;
    ComplexMatrix V,D;
//    StaticComplexMatrix.getEigenDecompose(m,&V,&D);
//    StaticComplexMatrix.print(V)END;
//    StaticComplexMatrix.print(D)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(StaticComplexMatrix.dot(V,D),StaticComplexMatrix.inverse(V)))END;
    //TODO not working properly for complex matrices
//    StaticComplexMatrix.getEigenDecomposeQR(m,&V,&D);
//    StaticComplexMatrix.print(V)END;
//    StaticComplexMatrix.print(D)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(StaticComplexMatrix.dot(V,D),StaticComplexMatrix.inverse(V)))END;
    ComplexMatrix U;
//    StaticComplexMatrix.getSVDDecompose(m,&U,&D,&V);
//    StaticComplexMatrix.print(U)END;
//    StaticComplexMatrix.print(D)END;
//    StaticComplexMatrix.print(V)END;
//    StaticComplexMatrix.print(StaticComplexMatrix.dot(StaticComplexMatrix.dot(U,D),StaticComplexMatrix.conjugate(V)))END;
    //TODO not working properly for rectangular matrices
    StaticComplexMatrix.getSVDDecomposeQR(m,&U,&D,&V);
    StaticComplexMatrix.print(U)END;
    StaticComplexMatrix.print(D)END;
    StaticComplexMatrix.print(V)END;
    StaticComplexMatrix.print(StaticComplexMatrix.dot(StaticComplexMatrix.dot(U,D),StaticComplexMatrix.conjugate(V)))END;



}

void linearTest(){
    float x[2][2]={
            {5,9},
            {5,9}
    };
    ComplexMatrix A = StaticComplexMatrix.getFromReal(2,2,x);
    StaticComplexMatrix.print(A)END;
    ComplexMatrix X = StaticLinear.solveHomogeneous(A);
    StaticComplexMatrix.print(X);
}

int main() {
//    complexTest();
//    polynomialTest();
//    matrixTest();
//    linearTest();
}
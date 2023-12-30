//
// Created by peter on 12/30/2023.
//

#ifndef C_LIB_POLYNOMIAL_H
#define C_LIB_POLYNOMIAL_H

#include "complex.h"

#define POLYNOMIAL_ITERATION 100

struct PolynomialControl {
    /**
     * This will calculate all the roots
     * @param n         : Degree of polynomial
     * @param function  : polynomial function
     * @param args      : Extra arguments may be needed
     * @param roots     : pointer to roots (array of complex number of size n)
     */
    void (*getRoots)(int n, Complex (*function)(Complex x, void **args), void **args, Complex *roots);
};

extern struct PolynomialControl StaticPolynomial;


#endif //C_LIB_POLYNOMIAL_H

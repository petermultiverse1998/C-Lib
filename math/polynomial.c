//
// Created by peter on 12/30/2023.
//

#include "polynomial.h"
#include "math.h"

#define c StaticComplex

/**
 * This will calculate all the roots
 * @param n         : Degree of polynomial
 * @param function  : polynomial function
 * @param args      : Extra arguments may be needed
 * @param roots     : pointer to roots (array of complex number of size n)
 */
static void getRoots(int n,Complex (*function)(Complex x,void** args),void** args,Complex* roots){
    Complex *z=roots;
    for (int k = 0; k < n; k++)
        z[k] = c.getFromPolar(1.0f,2.0f*(float)M_PI*(float)k/(float)n);
    for (int iteration = 0; iteration < 10; iteration++) {
        for (int k = 0; k < n; k++) {
            Complex x = z[k];
            Complex product = c.get(1,0);
            for (int j = 0; j <n ; j++) {
                if(j==k)
                    continue;
                product = c.dot(product,c.sub(x,z[j]));
            }

            z[k] = function(x,args);
            z[k] = c.div(function(x,args),product);
            z[k] = c.sub(x,z[k]);
        }
    }
}


struct PolynomialControl StaticPolynomial = {
        .getRoots = getRoots
};
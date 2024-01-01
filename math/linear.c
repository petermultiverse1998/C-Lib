//
// Created by peter on 12/30/2023.
//

#include "linear.h"

#define m StaticComplexMatrix
#define c StaticComplex

/**
 * This will solve homogeneous linear equation (AX=0)
 * @param augmentedMatrix   : Augmented matrix (A)
 * @return                  : Solved variables (X)
 */
static ComplexMatrix solveHomogeneous(ComplexMatrix augmentedMatrix){
    ComplexMatrix X = m.null(augmentedMatrix.row,1);
    ComplexMatrix Q,R;
    m.getQRWithGiven(augmentedMatrix,&Q,&R);

    int n = X.row;
    if(c.mag(R.a[n-1][n-1])>COMPLEX_EPSILON)
        X.a[n-1][0] = c.get(1,0);
    for (int k = n-2; k >=0; k--) {
        Complex sum = c.getFromReal(0);
        for (int j = k+1; j < n; j++)
            sum=c.sub(sum,c.dot(R.a[k][j],X.a[j][0]));
        if(c.mag(R.a[k][k])<COMPLEX_EPSILON)
            X.a[k][0]=c.get(1.0f,0.0f);
        else
            X.a[k][0]=c.div(sum,R.a[k][k]);
    }
    return X;
}

struct LinearControl StaticLinear = {
        .solveHomogeneous = solveHomogeneous
};
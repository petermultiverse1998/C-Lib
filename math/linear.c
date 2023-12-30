//
// Created by peter on 12/30/2023.
//

#include <stdio.h>
#include "linear.h"

/**
 * This will solve homogeneous linear equation (AX=0)
 * @param augmentedMatrix   : Augmented matrix (A)
 * @return                  : Solved variables (X)
 */
static ComplexMatrix solveHomogeneous(ComplexMatrix augmentedMatrix){
    ComplexMatrix X = StaticComplexMatrix.null(augmentedMatrix.row,1);
    ComplexMatrix Q,R;
    StaticComplexMatrix.getQRWithGiven(augmentedMatrix,&Q,&R);

    int n = X.row;
    if(StaticComplex.mag(R.a[n-1][n-1])>COMPLEX_EPSILON)
        X.a[n-1][0] = StaticComplex.get(1,0);
    for (int k = n-2; k >=0; k--) {
        Complex sum = StaticComplex.getFromReal(0);
        for (int j = k+1; j < n; j++)
            sum=StaticComplex.sub(sum,StaticComplex.dot(R.a[k][j],X.a[j][0]));
        if(StaticComplex.mag(R.a[k][k])<COMPLEX_EPSILON)
            X.a[k][0]=StaticComplex.get(1.0f,0.0f);
        else
            X.a[k][0]=StaticComplex.div(sum,R.a[k][k]);
    }
    return X;
}

struct LinearControl StaticLinear = {
        .solveHomogeneous = solveHomogeneous
};
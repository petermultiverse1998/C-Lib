//
// Created by peter on 12/30/2023.
//

#ifndef C_LIB_LINEAR_H
#define C_LIB_LINEAR_H

#include "complex_matrix.h"

struct LinearControl {
    /**
     * This will solve homogeneous linear equation (AX=0)
     * @param augmentedMatrix   : Augmented matrix (A)
     * @return                  : Solved variables (X)
     */
    ComplexMatrix (*solveHomogeneous)(ComplexMatrix augmentedMatrix);
};

extern struct LinearControl StaticLinear;


#endif //C_LIB_LINEAR_H

//
// Created by peter on 6/23/2023.
//

#ifndef C_LIB_HPF_H
#define C_LIB_HPF_H

typedef struct{
    float dt;
    int size;

    float tau;

    float *x;
    float *y;
}HPF;

struct HPFControl{
    /**
     * This will create object of HPF and return its pointer (It allocate memory so should be free using free)
     * @param dt    : Sampling interval
     * @param size  : size of inputs or outputs
     * @param tau   : Time constant or inverse of cutoff frequency
     * @return      : HPF pointer
     */
    HPF* (*new)(float dt,int size,float tau);

    /**
     * This will process and should be called in each step
     * @param hpf   : HPF
     * @param x     : inputs or errors
     * @return      : Same HPF
     */
    HPF* (*process)(HPF* hpf, const float* x);

    /**
     * This print the HPF
     * @param hpf   : HPF
     */
    void (*print)(HPF* hpf);

    /**
     * It free the allocated memory of HPF
     * @param hpfPtr    : pointer to HPF
     */
    void (*free)(HPF** hpfPtr);

    /**
     * This return allocated memory till now
     * @return  : Allocated memories
     */
    int (*getAllocatedMemories)();
};

extern struct HPFControl StaticHPF;

#endif //C_LIB_HPF_H

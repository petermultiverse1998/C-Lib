//
// Created by peter on 6/23/2023.
//

#ifndef C_LIB_LPF_H
#define C_LIB_LPF_H

typedef struct{
    float dt;
    int size;

    float tau;

    float *x;
    float *y;
}LPF;

struct LPFControl{
    /**
     * This will create object of LPF and return its pointer (It allocate memory so should be free using free)
     * @param dt    : Sampling interval
     * @param size  : size of inputs or outputs
     * @param tau   : Time constant or inverse of cutoff frequency
     * @return      : LPF pointer
     */
    LPF* (*new)(float dt,int size,float tau);

    /**
     * This will process and should be called in each step
     * @param lpf   : LPF
     * @param x     : inputs or errors
     * @return      : Same LPF
     */
    LPF* (*process)(LPF* lpf, const float* x);

    /**
     * This print the LPF
     * @param lpf   : LPF
     */
    void (*print)(LPF* lpf);

    /**
     * It free the allocated memory of LPF
     * @param lpfPtr    : pointer to LPF
     */
    void (*free)(LPF** lpfPtr);

    /**
     * This return allocated memory till now
     * @return  : Allocated memories
     */
    int (*getAllocatedMemories)();
};

extern struct LPFControl StaticLPF;



#endif //C_LIB_LPF_H

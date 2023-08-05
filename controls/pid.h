//
// Created by peter on 6/22/2023.
//

#ifndef C_LIB_PID_H
#define C_LIB_PID_H

typedef struct{
    float dt;
    int size;

    float Kp;
    float Ki;
    float Kd;
    float min;
    float max;

    float *x_prev;
    float *x;
    float *y;
}PID;

struct PIDControl{
    /**
     * This will create object of PID and return its pointer (It allocate memory so should be free using free)
     * @param dt    : Sampling interval
     * @param size  : size of inputs or outputs
     * @param Kp    : Proportional Constant
     * @param Ki    : Integrator Constant
     * @param Kd    : Derivative Constant
     * @param min   : Minimum Output
     * @param max   : Maximum Output
     * @return      : PID pointer
     */
    PID* (*new)(float dt,int size,float Kp,float Ki,float Kd,float min,float max);

    /**
     * This will process and should be called in each step
     * @param pid   : PID
     * @param x     : inputs or errors
     * @return      : Same PID
     */
    PID* (*process)(PID* pid, const float* x);

    /**
     * This print the PID
     * @param pid   : PID
     */
    void (*print)(PID* pid);

    /**
     * It free the allocated memory of PID
     * @param pidPtr    : pointer to PID
     */
    void (*free)(PID** pidPtr);

    /**
     * This return allocated memory till now
     * @return  : Allocated memories
     */
    int (*getAllocatedMemories)();
};

extern struct PIDControl StaticPID;

#endif //C_LIB_PID_H

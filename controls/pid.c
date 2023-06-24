//
// Created by peter on 6/22/2023.
//

#include <malloc.h>
#include <stdio.h>
#include "pid.h"

static int allocatedMemory = 0;

/**
 * It allocates the memory and return pointer to it
 * @param sizeInByte    : Size in bytes
 * @return              : Pointer to allocated memory
 *                      : NULL if there exist no memory for allocation
 */
static void *allocateMemory(int sizeInByte) {
    void* ptr = malloc(sizeInByte);
    if(ptr!=NULL)
        allocatedMemory+=sizeInByte;
    return ptr;
}

/**
 * It free the allocated memory
 * @param pointer       : Pointer to allocated Memory
 * @param sizeInByte    : Size to be freed
 * @return              : 1 for success (OR) 0 for failed
 */
static int freeMemory(void *pointer, int sizeInByte) {
    free(pointer);
    allocatedMemory-=sizeInByte;
    return 1;
}

/**
 * This will create object of PID and return its pointer (It allocate memory so should be free using free)
 * @param dt    : Sampling interval
 * @param size  : size of inputs or outputs
 * @param Kp    : Proportional Constant
 * @param Ki    : Integrator Constant
 * @param Kd    : Derivative Constant
 * @return      : PID pointer
 */
static PID* new(float dt,int size,float Kp,float Ki,float Kd){
    float a = Kp+Ki*dt+Kd/dt;
    float b = -Kp-2*Kd/dt;
    float c = Kd/dt;

    PID *pid = allocateMemory(sizeof(PID));

    pid->dt = dt;
    pid->size = size;

    pid->a = a;
    pid->b = b;
    pid->c = c;

    pid->x_prev = allocateMemory(size);
    pid->x = allocateMemory(size);
    pid->y = allocateMemory(size);
    for (int i = 0; i < size; ++i) {
        pid->x_prev[i]=0;
        pid->x[i]=0;
        pid->y[i]=0;
    }

    return pid;
}

/**
 * This will process and should be called in each step
 * @param pid   : PID
 * @param x     : inputs or errors
 * @return      : Same PID
 */
static PID* process(PID* pid, const float* xs){
    if(pid==NULL)
        return NULL;
    if(xs==NULL)
        return NULL;

    float y,x,x_prev;
    int size = pid->size;

    float a = pid->a;
    float b = pid->b;
    float c = pid->c;
    for (int i = 0; i < size; ++i) {
        y = pid->y[i];
        x = pid->x[i];
        x_prev = pid->x_prev[i];

        y = y + a*xs[i]+b*x+c*x_prev;

        pid->x_prev[i] = pid->x[i];
        pid->x[i] = xs[i];
        pid->y[i] = y;
    }
    return pid;
}

/**
 * This print the PID
 * @param pid   : PID
 */
static void print(PID* pid){
    if(pid==NULL)
        printf("PID: NULL\n");
    float Kd = pid->dt*pid->c;
    float Kp = -(pid->b+2*Kd/pid->dt);
    float Ki = (pid->a-Kp-Kd/pid->dt)/pid->dt;

    printf("PID: dt = %f, size = %d, Kp = %f, Ki = %f, Kd = %f\n",pid->dt,pid->size,Kp,Ki,Kd);
}

/**
 * It free the allocated memory of PID
 * @param pidPtr    : pointer to PID
 */
static void freePID(PID** pidPtr){
    PID* pid = *pidPtr;
    if(pid==NULL)
        return;
    if(pid->x_prev!=NULL)
        freeMemory(pid->x_prev, pid->size);
    if(pid->x!=NULL)
        freeMemory(pid->x, pid->size);
    if(pid->y!=NULL)
        freeMemory(pid->y, pid->size);
    freeMemory(pid, sizeof(PID));
    *pidPtr=NULL;
}

/**
 * This return allocated memory till now
 * @return  : Allocated memories
 */
static int getAllocatedMemories(){
    return allocatedMemory;
}

struct PIDControl StaticPID = {
        .new = new,
        .process = process,
        .print = print,
        .free = freePID,
        .getAllocatedMemories = getAllocatedMemories
};

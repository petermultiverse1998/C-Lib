//
// Created by peter on 6/23/2023.
//

#include "lpf.h"
#include <malloc.h>
#include <stdio.h>
#include <math.h>

static int allocatedMemory = 0;

/**
 * It allocates the memory and return pointer to it
 * @param sizeInByte    : Size in bytes
 * @return              : Pointer to allocated memory
 *                      : NULL if there exist no memory for allocation
 */
static void *allocateMemory(int sizeInByte) {
    void *ptr = malloc(sizeInByte);
    if (ptr != NULL)
        allocatedMemory += sizeInByte;
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
    allocatedMemory -= sizeInByte;
    return 1;
}

/**
 * This will create object of LPF and return its pointer (It allocate memory so should be free using free)
 * @param dt    : Sampling interval
 * @param size  : size of inputs or outputs
 * @param tau   : Time constant or inverse of cutoff frequency
 * @return      : LPF pointer
 */
static LPF *new(float dt, int size, float tau) {
    LPF *lpf = allocateMemory(sizeof(LPF));

    lpf->dt = dt;
    lpf->size = size;
    lpf->tau = tau;

    lpf->x = allocateMemory(size);
    lpf->y = allocateMemory(size);
    for (int i = 0; i < size; ++i) {
        lpf->x[i]=0;
        lpf->y[i]=0;
    }

    return lpf;
}

/**
 * This will process and should be called in each step
 * @param lpf   : LPF
 * @param xs    : inputs
 * @return      : Same LPF
 */
static LPF *process(LPF *lpf, const float *xs) {
    if (lpf == NULL)
        return NULL;
    if (xs == NULL)
        return NULL;

    float y;
    float dt = lpf->dt;
    int size = lpf->size;

    float tau = lpf->tau;
    for (int i = 0; i < size; ++i) {
        y = lpf->y[i];

        y = (y + dt * xs[i] / tau) / (1 + dt / tau);//Difference form
//        y = expf(-dt/tau)*(y-xs[i])+xs[i];//Exponential form

        lpf->y[i] = y;
        lpf->x[i] = xs[i];
    }
    return lpf;
}

/**
 * This print the LPF
 * @param lpf   : LPF
 */
static void print(LPF *lpf) {
    if (lpf == NULL)
        printf("LPF: NULL\n");
    printf("LPF: dt = %f, size = %d, tau = %f\n", lpf->dt, lpf->size, lpf->tau);
}

/**
 * It free the allocated memory of LPF
 * @param lpfPtr    : pointer to LPF
 */
static void freeLPF(LPF **lpfPtr) {
    LPF *lpf = *lpfPtr;
    if (lpf == NULL)
        return;
    if (lpf->x != NULL)
        freeMemory(lpf->x, lpf->size);
    if (lpf->y != NULL)
        freeMemory(lpf->y, lpf->size);
    freeMemory(lpf, sizeof(LPF));
    *lpfPtr = NULL;
}

/**
 * This return allocated memory till now
 * @return  : Allocated memories
 */
static int getAllocatedMemories() {
    return allocatedMemory;
}

struct LPFControl StaticLPF = {
        .new = new,
        .process = process,
        .print = print,
        .free = freeLPF,
        .getAllocatedMemories = getAllocatedMemories
};

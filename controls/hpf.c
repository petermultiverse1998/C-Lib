//
// Created by peter on 6/23/2023.
//

#include "hpf.h"
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
 * This will create object of HPF and return its pointer (It allocate memory so should be free using free)
 * @param dt    : Sampling interval
 * @param size  : size of inputs or outputs
 * @param tau   : Time constant or inverse of cutoff frequency
 * @return      : HPF pointer
 */
static HPF *new(float dt, int size, float tau) {
    HPF *hpf = allocateMemory(sizeof(HPF));

    hpf->dt = dt;
    hpf->size = size;
    hpf->tau = tau;

    hpf->x = allocateMemory(size);
    hpf->y = allocateMemory(size);
    for (int i = 0; i < size; ++i) {
        hpf->x[i]=0;
        hpf->y[i]=0;
    }

    return hpf;
}

/**
 * This will process and should be called in each step
 * @param hpf   : HPF
 * @param xs    : inputs
 * @return      : Same HPF
 */
static HPF *process(HPF *hpf, const float *xs) {
    if (hpf == NULL)
        return NULL;
    if (xs == NULL)
        return NULL;

    float y,x,q;
    float dt = hpf->dt;
    int size = hpf->size;

    float tau = hpf->tau;
    for (int i = 0; i < size; ++i) {
        y = hpf->y[i];
        x = hpf->x[i];

        y = (y + xs[i] - x) / (1 + dt / tau);//Difference form

        //Exponential form
//        q = (xs[i]-x)/dt;
//        y = expf(-dt/tau)*(y-q)+q;

        hpf->y[i] = y;
        hpf->x[i] = xs[i];
    }
    return hpf;
}

/**
 * This print the HPF
 * @param hpf   : HPF
 */
static void print(HPF *hpf) {
    if (hpf == NULL)
        printf("HPF: NULL\n");
    printf("HPF: dt = %f, size = %d, tau = %f\n", hpf->dt, hpf->size, hpf->tau);
}

/**
 * It free the allocated memory of HPF
 * @param hpfPtr    : pointer to HPF
 */
static void freeHPF(HPF **hpfPtr) {
    HPF *hpf = *hpfPtr;
    if (hpf == NULL)
        return;
    if (hpf->x != NULL)
        freeMemory(hpf->x, hpf->size);
    if (hpf->y != NULL)
        freeMemory(hpf->y, hpf->size);
    freeMemory(hpf, sizeof(HPF));
    *hpfPtr = NULL;
}

/**
 * This return allocated memory till now
 * @return  : Allocated memories
 */
static int getAllocatedMemories() {
    return allocatedMemory;
}

struct HPFControl StaticHPF = {
        .new = new,
        .process = process,
        .print = print,
        .free = freeHPF,
        .getAllocatedMemories = getAllocatedMemories
};

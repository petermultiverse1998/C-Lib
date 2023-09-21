//
// Created by peter on 9/20/2023.
//

#include "fft.h"
#include "math.h"
#include "stdio.h"

/**
 * This gives the Discrete Fourier Transform
 * @param N         : Sample Size
 * @param x         : Samples
 * @param dft_m     : Magnitudes
 * @param dft_ph    : Phases
 */
static void dft(int N,const float x[],float dft_m[],float dft_ph[]){
    // X_k = sum_n^(N-1) x_n*e^(-2*pi*i*k*n/N)
    for (int k = 0; k < N; ++k) {
        float real = 0;
        float img = 0;
        for (int n = 0; n < N; ++n) {
            real+=x[n]* cosf(2.0f*(float)M_PI*(float)k*(float)n/(float)N);
            img-=x[n]* sinf(2.0f*(float)M_PI*(float)k*(float)n/(float)N);
        }
        if(dft_m!=NULL)
            dft_m[k]= sqrtf(real*real+img*img);
        if(dft_ph!=NULL)
            dft_ph[k] = atan2f(img,real);
    }
}

/**
 * This gives the Inverse Discrete Fourier Transform
 * @param N         : Sample Size
 * @param dft_m     : Magnitudes
 * @param dft_ph    : Phases
 * @param x         : Constructed signals
 */
static void idft(int N,const float dft_m[],const float dft_ph[], float x[]){
    // x_k = (1/N)sum_n^(N-1) X_n*e^(2*pi*i*k*n/N)
    for (int k = 0; k < N; ++k) {
        float real = 0;
//        float img = 0;
        for (int n = 0; n < N; ++n) {
            real+=dft_m[n] * cosf(dft_ph[n]+2.0f*(float)M_PI*(float)k*(float)n/(float)N);
//            img+=dft_m[n]* sinf(dft_ph[n]+2.0f*(float)M_PI*(float)k*(float)n/(float)N);
        }
        real/=(float )N;
//        img/=(float )N;
        x[k]=real;
    }
}

typedef struct {
    float real;
    float img;
} Complex;

static void fftProcess(int N, const float x[], int origin, int spacing, Complex X[]) {
    if (N == 1) {
        X[0].real = x[origin];
        X[0].img = 0;
        return;
    }
    Complex X_even[N/2];
    Complex X_odd[N/2];

    fftProcess(N / 2, x, origin, 2 * spacing, X_even);
    fftProcess(N / 2, x, spacing + origin, 2 * spacing, X_odd);

    for (int k = 0; k < N / 2; k++) {
        Complex p = X_even[k];
        Complex w_N_k = {cosf(2.0f * (float) M_PI * (float) k / (float) N),
                         -sinf(2.0f * (float) M_PI * (float) k / (float) N)};
        Complex q = {w_N_k.real * X_odd[k].real - w_N_k.img * X_odd[k].img,
                     w_N_k.real * X_odd[k].img + w_N_k.img * X_odd[k].real};

        X[k].real = p.real + q.real;
        X[k].img = p.img + q.img;

        X[k + N / 2].real = p.real - q.real;
        X[k + N / 2].img = p.img - q.img;
    }
}

/**
 * This gives the Discrete Fourier Transform using FFT
 * @param N         : Sample Size
 * @param x         : Samples
 * @param dft_m     : Magnitudes
 * @param dft_ph    : Phases
 */
static void fft(int N,const float x[],float dft_m[],float dft_ph[]){
    Complex C[N];
    fftProcess(N,x,0,1,C);
    for(int i=0;i<N;i++) {
        if(dft_m!=NULL)
            dft_m[i] = sqrtf(C[i].real * C[i].real + C[i].img * C[i].img);
        if(dft_ph!=NULL)
            dft_ph[i] = atan2f(C[i].img,C[i].real);
    }
}

static void ifftProcess(int N,const Complex X[], int origin, int spacing,Complex x[]) {
    if (N == 1) {
        x[0].real = X[origin].real;
        x[0].img = X[origin].img;
        return;
    }
    Complex x_even[N/2];
    Complex x_odd[N/2];

    ifftProcess(N / 2, X,origin, 2 * spacing, x_even);
    ifftProcess(N / 2 , X,spacing + origin, 2 * spacing, x_odd);

    for (int k = 0; k < N / 2; k++) {
        Complex p = {x_even[k].real/2.0f,x_even[k].img/2.0f};
        Complex w_N_k = {cosf(2.0f * (float) M_PI * (float) k / (float) N),
                         sinf(2.0f * (float) M_PI * (float) k / (float) N)};
        Complex q = {(w_N_k.real * x_odd[k].real - w_N_k.img * x_odd[k].img)/2.0f,
                     (w_N_k.real * x_odd[k].img + w_N_k.img * x_odd[k].real)/2.0f};

        x[k].real = p.real + q.real;
        x[k].img = p.img + q.img;

        x[k + N / 2].real = p.real - q.real;
        x[k + N / 2].img = p.img - q.img;
    }
}

/**
 * This gives the inverse Discrete Fourier Transform using FFT
 * @param N         : Sample Size
 * @param dft_m     : Magnitudes
 * @param dft_ph    : Phases
 * @param x         : Samples
 */
static void ifft(int N,const float dft_m[],float dft_ph[],float x[]){
    Complex X[N];
    Complex C[N];
    for (int i = 0; i < N; ++i) {
        X[i].real = dft_m[i]* cosf(dft_ph[i]);
        X[i].img = dft_m[i]* sinf(dft_ph[i]);
    }

    ifftProcess(N,X,0,1,C);
    for(int i=0;i<N;i++)
        x[i]=C[i].real;
}


struct FFTControl StaticFFT = {
        .dft = dft,
        .idft = idft,
        .fft = fft,
        .ifft = ifft
};
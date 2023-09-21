//
// Created by peter on 9/20/2023.
//

#ifndef C_LIB_FT_H
#define C_LIB_FT_H

struct FFTControl{
    /**
     * This gives the Discrete Fourier Transform
     * @param N         : Sample Size
     * @param x         : Samples
     * @param dft_m     : Magnitudes
     * @param dft_ph    : Phases
     */
    void (*dft)(int N,const float x[],float dft_m[],float dft_ph[]);

    /**
     * This gives the Inverse Discrete Fourier Transform
     * @param N         : Sample Size
     * @param dft_m     : Magnitudes
     * @param dft_ph    : Phases
     * @param x         : Constructed signals
     */
    void (*idft)(int N,const float dft_m[],const float dft_ph[], float x[]);

    /**
     * This gives the Discrete Fourier Transform using FFT
     * @param N         : Sample Size
     * @param x         : Samples
     * @param dft_m     : Magnitudes
     * @param dft_ph    : Phases
     */
    void (*fft)(int N,const float x[],float dft_m[],float dft_ph[]);

    /**
     * This gives the inverse Discrete Fourier Transform using FFT
     * @param N         : Sample Size
     * @param dft_m     : Magnitudes
     * @param dft_ph    : Phases
     * @param x         : Samples
     */
    void (*ifft)(int N,const float dft_m[],float dft_ph[],float x[]);
};

extern struct FFTControl StaticFFT;
#endif //C_LIB_FT_H

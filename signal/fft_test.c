//
// Created by peter on 9/20/2023.
//
#include "fft.h"
#include "stdio.h"

void print(int N,const float array[]){
    printf("=>");
    for (int i = 0; i < N; ++i)
        printf("%.2f ",array[i]);
    printf("\n");
}

int main(){
    int N = 4;
    float x[]={1,2,3,4,5,6,7,8};
    float X[N];
    float ph[N];
    float y[N];

    print(N,x);

//    StaticFFT.dft(N,x,X,ph);
//    printf("DFT\n");
//    print(N,X);
//    print(N,ph);
//    printf("IDFT\n");
//    StaticFFT.idft(N,X,ph,y);
//    print(N,y);

    StaticFFT.fft(N,x,X,ph);
    printf("FFT\n");
    print(N,X);
    print(N,ph);
    printf("IFFT\n");
    StaticFFT.ifft(N,X,ph,y);
    print(N,y);


}
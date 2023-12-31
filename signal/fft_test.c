//
// Created by peter on 9/20/2023.
//
#include <math.h>
#include "fft.h"
#include "stdio.h"
#include "windows.h"



static void copyToClipboard(const char *output) {
//    const char* output = "Test\0abc";
    const size_t len = strlen(output) + 1;
    HGLOBAL hMem = GlobalAlloc(GMEM_MOVEABLE, len);
    memcpy(GlobalLock(hMem), output, len);
    GlobalUnlock(hMem);
    OpenClipboard(0);
    EmptyClipboard();
    SetClipboardData(CF_TEXT, hMem);
    CloseClipboard();
}

//float x[]={5,2,4,0,4,3,1,2,1,7,3,2,4,7,2,0,4,0,0,1,1,7,2,0,9,3,3,7,3,2,0,4,2,0,0,1,1,0,4,3,2,0,3,0,2,2,0,1,3,2,0,1,2,0,1,3,2,2,0,1,4,1,0,0,2,3,0,3,0,0,2,2,0,2,1,1,2,2,0,0,0,1,0,3,1,1,2,1,0,0,0,0,1,2,4,0,0,0,2,2,0,2,0,0,1,0,1,0,1,3,2,1,1,2,0,0,4,0,8,2,1,0,2,0,2,1,0,0,2,0,1,2,6,0,2,0,1,1,2,1,0,1,1,0,1,0,0,2,1,2,0,2,0,0,0,3,8,1,0,1,0,0,2,0,14,1,0,1,0,1,0,0,4,1,0,1,1,0,0,0,3,0,2,2,1,0,0,2,1,3,6,3,10,7,3,5,5,8,6,8,10,12,13,13,16,15,12,18,22,22,22,26,27,30,28,27,27,32,30,30,33,32,32,29,34,34,36,33,36,34,33,35,33,32,33,34,36,35,32,34,35,33,36,35,35,36,33,33,39,35,32,28,36,35,40,36,34,33,35,32,33,36,34,35,34,35,35,34,30,33,33,40,34,34,34,34,37,35,33,32,33,34,34,37,34,33,29,33,37,33,36,36,34,36,33,33,24,34,35,33,37,34,35,34,34,32,34,39,36,34,37,36,34,33,34,34,35,34,33,39,34,30,36,33,34,39,36,36,32,36,34,33,31,32,35,34,32,34,34,35,32,30,34,35,34,35,34,32,33,34,34,36,34,36,34,35,34,34,32,40,34,34,31,35,34,35,39,35,34,33,32,36,33,32,37,34,33,35,32,36,34,32,35,34,34,33,29,36,34,39,34,35,34,33,34,34,35,42,36,33,36,34,34,30,36,34,34,39,35,33,33,33,33,35,34,33,35,35,34,35,34,33,33,36,35,36,38,34,33,33,33,32,35,36,35,37,33,34,37,28,34,32,30,24,21,17,14,15,15,13,14,13,10,10,9,10,5,6,8,8,9,7,6,8,5,6,4,8,3,5,4,5,3,4,3,3,2,10,5,2,3,3,3,2,0,4,4,8,5,2,2,3,2,1,3,1,9,3,2,1,4,0,1,3,2,4,1,2,2,7,3,1,2,0,4,2,1,2,2,0,2,4,0,3,1,1,2,2,1,0,1,3,0,3,2,1,0,2,0,2,3,2,0,1,2,0,2,2,1,0,4,1,4,2,2,0,1,2,0,4,0,0,4,2,0,0,4,0,3,3,1,0,3,0,0,2,3,7,0,3,1,0,1,1,2,1,5,0,1,0,2,0,0,3,1,0,1,1,1,0,2,2,0,8,2,2,0,1,2,0,0,1,0,1,3,0,0,3,0,0,5,2,0,4,0,0,1,0,0,3,6,1,5,1,0,0,1,4,1,2,0,2,4,6,0,0,1,2,0,0,2,0,1,2,1,0,0,1,3,2,1,1,1,2,0,1,7,0,2,0,1,1,1,1,3,5,1,1,1,0,0,0,1,0,3,1,1,0,1,0,1,1,2,0,8,2,0,0,2,1,0,0,2,2,0,0,3,0,0,0,3,2,5,2,3,2,0,5,3,6,15,12,7,5,12,9,10,12,14,15,14,20,18,19,20,21,30,27,26,27,29,27,28,32,29,30,33,30,32,35,31,33,32,30,32,33,33,39,34,32,32,35,28,36,34,37,35,34,34,33,33,34,35,33,34,34,34,33,35,34,33,34,36,34,34,34,35,30,36,33,33,31,35,34,30,36,33,32,34,34,34,33,30,34,33,34,34,27,32,36,35,34,34,30,31,34,31,34,34,34,36,34,34,29,34,33,40,35,34,35,34,33,32,36,32,33,33,35,34,34,36,34,33,33,35,33,36,34,34,35,33,33,28,33,35,32,35,35,36,33,35,34,36,37,34,33,34,34,34,33,35,33,35,34,34,34,35,33,33,34,34,36,34,34,33,34,34,32,34,28,36,36,31,33,34,33,32,28,35,35,40,35,34,33,32,32,34,32,37,35,33,34,34,34,32,35,33,34,32,35,35,30,34,33,33,35,34,33,36,33,33,30,35,34,34,35,37,36,34,34,32,32,33,34,34,38,36,33,33,33,34,33,34,39,34,33,33,36,34,33,36,33,35,33,32,27,23,20,20,16,16,13,11,15,9,10,9,10,9,7,12,7,8,7,5,5,3,4,6,3,13,5,4,4,4,4,2,4,19,6,2,4,4,2,2,5,3,2,5,3,1,11,2,2,3,4,1,4,3,4,3,3,1,9};

int main(){
    int N = 1024;//Samples
    float x[N];
    float dt = 1.0f/(float)N;//sampling time
    for (int i = 0; i < N; ++i)
        x[i] = sinf(2.0f*(float)M_PI*dt*(float)i);//Sine Wave

    float X[N];
    float ph[N];

    StaticFFT.fft(N,x,X,ph);//FFT calculate
    for(int i=0;i<N/2;i++){
        float temp=X[i]/(float)N;
        if(temp<=0.4f) {
            X[i] = 0;
            X[N-i]=0;
        }
    }
    StaticFFT.ifft(N,X,ph,x);

    char str[10024];
    int ptr = 0;
    for (int i = 0; i < sizeof(X) / 4; ++i)
        ptr += sprintf(str + ptr, "%0.4f \n", x[i]);
//        ptr += sprintf(str + ptr, "%0.4f \n",X[i]);
//        ptr += sprintf(str + ptr, "%0.4f \n",ph[i]);
    copyToClipboard(str);
}
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

int main(){
    int N = 256;//Samples
    float x[N];
    float dt = 1.0f/(float)N;//sampling time
    for (int i = 0; i < N; ++i)
        x[i] = sinf(2.0f*(float)M_PI*dt*(float)i);//Sine Wave

    float X[N];
    float ph[N];

    StaticFFT.fft(N,x,X,ph);

    char str[10024];
    int ptr = 0;
    for (int i = 0; i < sizeof(X) / 4; ++i)
        ptr += sprintf(str + ptr, "%f \n", x[i]);
//        ptr += sprintf(str + ptr, "%0.4f \n",X[i]);
//        ptr += sprintf(str + ptr, "%0.4f \n",ph[i]);
    copyToClipboard(str);
}
//
// Created by peter on 6/23/2023.
//

#include "hpf.h"

#include <stdio.h>
#include <windows.h>
#include <math.h>
#include "../test/test.h"

static int testNew() {
    float dt = 0.01f;
    int size = 1;
    float tau = 1;
    HPF *hpf = StaticHPF.new(dt, size, tau);

    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != hpf);

    check = check && FLOAT_EQUALS(__FILE__, __LINE__, dt, hpf->dt);
    check = check && INT_EQUALS(__FILE__, __LINE__, size, hpf->size);
    check = check && FLOAT_EQUALS(__FILE__, __LINE__, tau, hpf->tau);
    check = check && BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != hpf->x);
    check = check && BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != hpf->y);

    printf("testNew\n");
    StaticHPF.free(&hpf);
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticHPF.getAllocatedMemories());
    return check;
}

static int testProcess() {
    float dt = 0.01f;
    int size = 1;
    float tau = 1;
    HPF *hpf = StaticHPF.new(dt, size, tau);

    float x[] = {1};
    float y[size];
    for (int i = 0; i < size; ++i) {
        /*
         * y = (y+x-x_prev)/(1+dt/T)
         */
        y[i] = 0;
        y[i] = (y[i] + x[i]-0) / (1 + dt / tau);
    }

    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != StaticHPF.process(hpf, x));
    check = check && FLOAT_ARRAY_EQUALS(__FILE__, __LINE__, y, hpf->y, size);

    printf("testProcess\n");
    StaticHPF.free(&hpf);
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticHPF.getAllocatedMemories());
    return check;
}

static int testFree() {
    float dt = 0.01f;
    int size = 1;
    float tau = 1;
    HPF *hpf = StaticHPF.new(dt, size, tau);

    StaticHPF.free(&hpf);
    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL == hpf);
    StaticHPF.free(&hpf);

    printf("testFree\n");
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticHPF.getAllocatedMemories());
    return check;
}

static void test() {
    Test test = StaticTest.new();
    StaticTest.addTask(&test, testNew);
    StaticTest.addTask(&test, testProcess);
    StaticTest.addTask(&test, testFree);

    StaticTest.run(&test);
}

static void copyToClipboard(const char* output){
//    const char* output = "Test\0abc";
    const size_t len = strlen(output) + 1;
    HGLOBAL hMem =  GlobalAlloc(GMEM_MOVEABLE, len);
    memcpy(GlobalLock(hMem), output, len);
    GlobalUnlock(hMem);
    OpenClipboard(0);
    EmptyClipboard();
    SetClipboardData(CF_TEXT, hMem);
    CloseClipboard();
}

static void demo() {
    float dt = 0.05f;
    float T = 10.0f;
    int N = (int)(T/dt);
    int size = 1;

    float x[N][size],y[N][size];
    for (int i = 0; i < N; ++i) {
//        x[i][0] = 1;//Step
//        x[i][0] = (i>=0 && i<=10)?1:0;//Pulse
//        x[i][0] = (i%20<10)?1:0;//Square Wave
//        x[i][0] = (i%20<10)?((float)(i%10)/10.0f):1-((float)(i%10)/10.0f);//Triangular Wave
//        x[i][0] = (float)(i%21)/20.0f;//Sawtooth Wave
        x[i][0] = sinf(0.5f*2.0f*(float)M_PI*dt*(float)i);//Sine Wave

    }

    HPF *hpf = StaticHPF.new(dt, size, 1);
    for (int i = 0; i < N; ++i) {
        StaticHPF.process(hpf,x[i]);
        y[i][0] = hpf->y[0];
    }

    char str[N*10];
    int ptr = 0;
    for (int i = 0; i < N; ++i)
        ptr+=sprintf(str+ptr, "%f \n",x[i][0]);
    copyToClipboard(str);

    StaticHPF.print(hpf);
    StaticHPF.free(&hpf);
}

int main() {
//    test();
    demo();

    return 0;
}

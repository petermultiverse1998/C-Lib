//
// Created by peter on 6/22/2023.
//
#include "pid.h"
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include "../test/test.h"

static int testNew() {
    float dt = 0.01f;
    int size = 1;

    float Kp = 1;
    float Ki = 1;
    float Kd = 1;

    float a = Kp+Ki*dt+Kd/dt;
    float b = -Kp-2*Kd/dt;
    float c = Kd/dt;


    PID *pid = StaticPID.new(dt, size, Kp,Ki,Kd);

    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != pid);

    check = check && FLOAT_EQUALS(__FILE__, __LINE__, dt, pid->dt);
    check = check && INT_EQUALS(__FILE__, __LINE__, size, pid->size);
    check = check && FLOAT_EQUALS(__FILE__, __LINE__, a, pid->a);
    check = check && FLOAT_EQUALS(__FILE__, __LINE__, b, pid->b);
    check = check && FLOAT_EQUALS(__FILE__, __LINE__, c, pid->c);
    check = check && BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != pid->x_prev);
    check = check && BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != pid->x);
    check = check && BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != pid->y);

    printf("testNew\n");
    StaticPID.free(&pid);
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticPID.getAllocatedMemories());
    return check;
}

static int testProcess() {
    float dt = 0.01f;
    int size = 1;

    float Kp = 1;
    float Ki = 1;
    float Kd = 1;

    float a = Kp+Ki*dt+Kd/dt;
    float b = -Kp-2*Kd/dt;
    float c = Kd/dt;


    PID *pid = StaticPID.new(dt, size, Kp,Ki,Kd);

    float x[] = {1};
    float y[size];
    for (int i = 0; i < size; ++i) {
        /*
         * y = y+ax+bx_prev+cx_prev_prev
         */
        y[i] = 0;
        y[i] = y[i]+a*x[i]+b*0+c*0;
    }

    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL != StaticPID.process(pid, x));
    check = check && FLOAT_ARRAY_EQUALS(__FILE__, __LINE__, y, pid->y, size);

    printf("testProcess\n");
    StaticPID.free(&pid);
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticPID.getAllocatedMemories());
    return check;
}

static int testFree() {
    float dt = 0.01f;
    int size = 1;

    float Kp = 1;
    float Ki = 1;
    float Kd = 1;

    PID *pid = StaticPID.new(dt, size, Kp,Ki,Kd);

    StaticPID.free(&pid);
    int check = BOOLEAN_IS_TRUE(__FILE__, __LINE__, NULL == pid);
    StaticPID.free(&pid);

    printf("testFree\n");
    check = check && INT_EQUALS(__FILE__, __LINE__, 0, StaticPID.getAllocatedMemories());
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


    PID *pid = StaticPID.new(dt, size, 0.001f,20,0);


    float y_ref[1],error[1];
    for (int i = 0; i < N; ++i) {
        y_ref[0] = x[i][0];
        error[0] = y_ref[0]-pid->y[0];
        StaticPID.process(pid,error);
        y[i][0] = pid->y[0];
    }

    char str[N*10];
    int ptr = 0;
    for (int i = 0; i < N; ++i)
        ptr+=sprintf(str+ptr, "%f \n",x[i][0]);
    copyToClipboard(str);

    StaticPID.print(pid);
    StaticPID.free(&pid);
}

int main() {
//    test();
    demo();

    return 0;
}

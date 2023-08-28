#include <stdio.h>
#include "test/test.h"

#define varToStr(var) #var

void synchronise(){
    //Make function busy
    static int isBusy = 0;
    while(isBusy);
    isBusy = 1;

    //To Work Here

    //Release busy flag
    isBusy = 0;
}


int main() {
    int x = 2;
//    RED;


    return 0;
}

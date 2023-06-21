#include <stdio.h>

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

    return 0;
}

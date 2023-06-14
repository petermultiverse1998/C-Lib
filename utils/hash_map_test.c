//
// Created by peter on 6/13/2023.
//

#include "hash_map.h"
#include <assert.h>
#include "sys/time.h"
#include "stdio.h"
#include "../test/test.h"

static long getTimeInMicros() {
    struct timeval val;
    gettimeofday(&val, NULL);
    return val.tv_usec;
}

static int test1(){
    HashMap *map = StaticHashMap.new();
    return BOOLEAN_IS_TRUE(__FILE__,__LINE__,NULL!=map);
}

void hash_test() {
    Test test = StaticTest.new();

    StaticTest.addTask(&test,test1);

    StaticTest.run(&test);
}
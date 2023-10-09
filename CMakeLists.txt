cmake_minimum_required(VERSION 3.25)
project(C_Lib C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
add_executable(
        ${PROJECT_NAME}
        test/test.c
        test/test.h

        utils/hash_map.c
        utils/hash_map.h
        utils/queue.c
        utils/queue.h
        utils/task.c
        utils/task.h
        utils/static_heap.c
        utils/static_heap.h
        utils/buddy_heap.c
        utils/buddy_heap.h

        controls/lpf.c
        controls/lpf.h
        controls/hpf.c
        controls/hpf.h
        controls/pid.c
        controls/pid.h

        signal/fft.c
        signal/fft.h

#        test/test_test.c
#
#        utils/hash_map_test.c
#        utils/queue_test.c
#        utils/task_test.c
#        utils/static_heap_test.c
        utils/buddy_heap_test.c
#
#        controls/lpf_test.c
#        controls/hpf_test.c
#        controls/pid_test.c
#
#        signal/fft_test.c

#        main.c
)
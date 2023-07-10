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

        controls/lpf.c
        controls/lpf.h
        controls/hpf.c
        controls/hpf.h
        controls/pid.c
        controls/pid.h

        #        test/test_test.c

        #        utils/hash_map_test.c
                utils/queue_test.c
        #        utils/task_test.c

        #        controls/lpf_test.c
        #        controls/hpf_test.c
        #        controls/pid_test.c

#        main.c
)
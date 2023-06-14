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

        utils/linked_list.c
        utils/linked_list.h

#        test/test_test.c
        utils/hash_map_test.c
#        main.c

)
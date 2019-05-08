cmake_minimum_required(VERSION 2.8)

project(superlu C)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../../Libs)

file(GLOB SUPERLU_SRC "*.c")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DNDEBUG -DPRNTlevel=0 -DAdd_")

if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -W0")
else()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
endif()

add_library(${PROJECT_NAME} STATIC ${SUPERLU_SRC})
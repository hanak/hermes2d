if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(lshape-square)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

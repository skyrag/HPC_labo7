cmake_minimum_required(VERSION 3.10)
project(FindKmers)

set(CMAKE_C_STANDARD 99)

add_executable(k-mer src/main.c)
add_executable(k-mer-parallel src/parallelmain.c)

# Compile the executable with optimizations by default
target_compile_options(k-mer PRIVATE -O0 -g)
target_compile_options(k-mer-parallel PRIVATE -O0 -g -fopenmp)


find_package(OpenMP REQUIRED)

target_link_libraries(k-mer-parallel PRIVATE OpenMP::OpenMP_C)
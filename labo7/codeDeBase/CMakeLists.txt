cmake_minimum_required(VERSION 3.10)
project(FindKmers)

set(CMAKE_C_STANDARD 99)

add_executable(k-mer src/main.c)

# Compile the executable with optimizations by default
target_compile_options(k-mer PRIVATE -O0 -g)

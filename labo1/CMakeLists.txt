cmake_minimum_required(VERSION 3.10)
project(ecg_dealination C)

set(CMAKE_C_STANDARD 99)

set(SOURCES
    src/main.c
    src/ecg_processing.c
    src/csv_reader.c
    src/json_writer.c
    src/ecg_utils.c
)

include_directories(include)

add_executable(ecg_dealination ${SOURCES})

add_executable(ecg_dealination_parallel src/parallelMain.c
    src/ecg_processing.c
    src/csv_reader.c
    src/json_writer.c
    src/ecg_utils.c
)

find_package(OpenMP REQUIRED)

target_link_libraries(ecg_dealination m)
target_link_libraries(ecg_dealination_parallel PRIVATE m OpenMP::OpenMP_C)

target_compile_options(ecg_dealination_parallel PRIVATE -O0 -fopenmp)
target_compile_options(ecg_dealination PRIVATE -O0)



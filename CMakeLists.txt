#cmake_minimum_required(VERSION 3.27)
cmake_minimum_required(VERSION 2.6)

#************************* GPU VERSION *******************************

#find_package(CUDA REQUIRED)
set(CUDA_DIR /usr/local/cuda)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CUDA_ARCHITECTURES native)
if(CMAKE_CUDA_COMPILER)

    project(difx_gpu C CXX CUDA)

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -g -Wall ${ADDITIONAL_C_FLAGS}")

    #find_package(CUDA 12.3 REQUIRED)
    find_package(CUDA REQUIRED)
    message(STATUS ">>> Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")

    set(
            CUDA_NVCC_FLAGS
            ${CUDA_NVCC_FLAGS};
            -std=c++11 -arch=sm_60 -Xptxas -O3 -restrict -D_FORCE_INLINES -D_DEBUG -lineinfo
    )

    file(
            GLOB
            source_files_GPU
            gkernel.cu
            main_gpu.cu
            openFiles.c
    )

    file(
            GLOB
            header_files
            gkernel.cuh
            openFiles.h
    )

    message("source_files_GPU : ${source_files_GPU}")
    #    message("header_files : ${header_files}")
    cuda_add_executable(difx_gpu ${source_files_GPU})
    target_link_libraries(difx_gpu cudart cusolver cufft m) # math, CUDA FFT and SVD

else()
    message(STATUS "No CUDA compiler found")
endif()

#************************* CPU VERSION *******************************

project(tests_difx C)

set(CMAKE_C_STANDARD 11)

set(IPP_DIR /opt/intel/oneapi/ipp/latest/lib/cmake/ipp)

FIND_PACKAGE(IPP REQUIRED)

if (IPP_FOUND)
add_definitions (-DHAVE_IPP)
set(IPP_INCLUDE_DIRS /opt/intel/oneapi/ipp/latest/include/)
set(IPP_LIBRARY_DIRS /opt/intel/oneapi/ipp/latest/lib/)

include_directories (${IPP_INCLUDE_DIRS})
include_directories(${IPP_INCLUDE_DIRS})
link_directories ( ${IPP_LIBRARY_DIRS} )

find_library(IPP_CORE libippcore.so PATHS /opt/intel/oneapi/ipp/latest/lib/intel64)

endif ()

add_executable(tests_difx main.c
        allocDataHost.h
        allocDataHost.c
        openFiles.h
        openFiles.c
        kernel.h
        kernel.c
        vectordefs.h
)

TARGET_LINK_LIBRARIES(tests_difx ${IPP_LIBRARIES} ${CMAKE_MATH_LIBRARIES} -lm)


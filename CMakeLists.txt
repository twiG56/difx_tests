cmake_minimum_required(VERSION 3.27)
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
        vectordefs.h)

TARGET_LINK_LIBRARIES(tests_difx ${IPP_LIBRARIES} ${CMAKE_MATH_LIBRARIES} -lm)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(g2o_slim_EDU)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
set(CMAKE_BUILD_TYPE Release)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#find macros and libraries
find_package(Eigen3 REQUIRED)
find_package(Boost  REQUIRED)
find_package(SuiteSparse REQUIRED)

#include directories
include_directories(
	src
	${EIGEN3_INCLUDE_DIR}
	${Boost_INCLUDE_DIR}
	${OpenCV_INCLUDE_DIRS}
	${CSPARSE_INCLUDE_DIR})


add_subdirectory(src)
add_subdirectory(unittests)


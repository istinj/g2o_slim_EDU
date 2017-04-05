CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(test_Cholesky)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_modules)
set(CMAKE_BUILD_TYPE Release)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

#find macros and libraries
find_package(Eigen3 REQUIRED
	Boost 1.54.0.1 REQUIRED
	G2O REQUIRED
	OpenCV REQUIRED
	CSparse REQUIRED)

#include directories
include_directories(
	include
	${EIGEN3_INCLUDE_DIR}
	${Boost_INCLUDE_DIR}
	${G2O_INCLUDE_DIR}
	${OpenCV_INCLUDE_DIRS}
	${CSPARSE_INCLUDE_DIR})

set(SRCS 
	src/test.cpp
	src/utilities.cpp
	src/SparseMatrix.cpp)

add_executable(test_Cholesky_app ${SRCS})
# Specify libraries to link a library or executable target against
target_link_libraries(test_Cholesky_app ${Boost_LIBRARIES})
set_target_properties( test_Cholesky_app PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )
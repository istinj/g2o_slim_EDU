project(tls)
cmake_minimum_required( VERSION 2.8 )
set(G2O_SLIM_EDU_ROOT /home/istin/source/g2o_slim_EDU)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} --std=gnu++11)
set(CMAKE_BUILD_TYPE Release)

#find macros and libraries
find_package(Eigen3 REQUIRED
	Boost 1.54.0.1 REQUIRED)

#include directories
include_directories(include
	${G2O_SLIM_EDU_ROOT}/sparse/include
	${EIGEN3_INCLUDE_DIR}
	${Boost_INCLUDE_DIR})
	
message(${G2O_SLIM_EDU_ROOT}/sparse/include)

#c++ executable
set(SRCS 
	src/tls_main.cpp
	src/SparseSolver.cpp
	src/Graph.cpp)
add_executable(tls_main_app ${SRCS})

# Specify libraries to link a library or executable target against
target_link_libraries(tls_main_app ${Boost_LIBRARIES})

set_target_properties( tls_main_app PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin )

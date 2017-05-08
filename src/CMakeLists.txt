add_library(g2o_slim_EDU_library SHARED
	defs.h
	edge.h							edge.cpp
	factor.h						factor.cpp
	graph.h							graph.cpp
	sparse_block_matrix.h			sparse_block_matrix.cpp
	sparse_optimizer.h				sparse_optimizer.cpp
	vertex.h						vertex.cpp
	workspace.h						workspace.cpp)
	
	
target_link_libraries(g2o_slim_EDU_library
	${Boost_LIBRARIES} ${CSPARSE_LIBRARY})

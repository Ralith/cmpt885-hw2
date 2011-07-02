add_executable(genetic-p
  genetic_tsp.cpp
  )

find_library(M_LIB m)
target_link_libraries(genetic-p ${M_LIB})

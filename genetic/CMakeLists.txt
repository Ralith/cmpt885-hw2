add_executable(genetic
  genetic_tsp.cpp
  )

find_library(M_LIB m)
target_link_libraries(genetic ${M_LIB})

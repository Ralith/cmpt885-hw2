add_executable(ant
  ant_colony_tsp.cpp
  )

find_library(M_LIB m)
target_link_libraries(ant ${M_LIB})

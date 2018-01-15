# everytime we have a new directory, we manually add the directory here
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/aqdata)

file(GLOB ALL_SRCS "*.cc")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

ADD_LIBRARY(_aqdata ${ALL_SRCS})

DEAL_II_SETUP_TARGET(_aqdata)

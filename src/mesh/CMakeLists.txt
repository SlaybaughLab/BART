# everytime we have a new directory, we manually add the 
# absolute directory here
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/mesh)

FILE(GLOB ALL_SRCS "*.cc")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# change the library names correspondingly for new folders
ADD_LIBRARY(_mesh ${ALL_SRCS})

DEAL_II_SETUP_TARGET(_mesh)

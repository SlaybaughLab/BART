#ifndef BART_SRC_UTILITY_BART_HASHER_H_
#define BART_SRC_UTILITY_BART_HASHER_H_

#include <utility>
#include <tuple>
#include <vector>

#include <boost/functional/hash.hpp>

//! Provides various utility functions for BART.
namespace butil {

//! Hashing struct for array<int> of any size, and vector<int>.
/*! Provides a hashing struct for array<int> of any size, and vector<int>;
 *  Suitable for use in an unordered_map. Uses the hashing functions of the boost
 *  library.
 *  Example:
 *  ~~~~~~~~~~{.cpp}
 *  #include "bart_hasher.h"
 *  #include <array>
 *  #include <iostream>
 *  #include <string>
 *  #include <unordered_map>
 *
 *  int main()
 *  {
 *    std::unordered_map<std::array<int, 2>, std::string, butil::Hasher>
 *      array_to_string_umap = {
 *                              {{0,0}, "all zeros"},
 *                              {{1,1}, "all ones"}
 *                             };
 *    std::cout << array_to_string_umap[{1,1}] << "\n"
 *              << array_to_string_umap[{0,0}] << std::endl;
 *    return 0;
 *  }
 *  ~~~~~~~~~~
 *  Output:
 *  ~~~~~~~~~~
 *  all ones
 *  all zeros
 *  ~~~~~~~~~~
 */    

struct Hasher {
  //! Converts an array<int, array_size> into a hashed value.
  template <std::size_t array_size>
  std::size_t operator()(const std::array<int, array_size>) const;
  
  //! Converts a vector<int> into a hashed value.
  std::size_t operator()(const std::vector<int>) const;
};

} //namespace butil

#endif

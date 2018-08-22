#ifndef BART_SRC_UTILITY_BART_HASHER_H_
#define BART_SRC_UTILITY_BART_HASHER_H_

#include <utility>
#include <tuple>
#include <vector>

#include <boost/functional/hash.hpp>

//! This namespace provides various utility functions for BART.
namespace butil {

//! Hashing struct for array<int> of any size, and vector<int>.
/*! Provides a hashing struct for array<int> of any size, and vector<int>;
 *  Suitable for use in an unordered_map. Uses the hashing functions of the boost
 *  library.
 *  Example:
 *    std::unordered_map<array<int, 2>, std::string, butil::Hasher> map =
 *      {{0, 0}, "all zeros"},
 *      {{1, 1}, "all ones"},};
 *    std::string value = map[{1,1}]; // Value is "all ones"
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

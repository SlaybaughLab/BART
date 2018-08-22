#include "bart_hasher.h"

namespace butil {

template <std::size_t array_size>
std::size_t Hasher::operator()(const std::array<int, array_size> array) const {
  std::size_t seed = 0;
  // Iterates over array values and combines their hashed values
  for (const int val : array)
    boost::hash_combine(seed, boost::hash_value(val));
  return seed;
}

std::size_t Hasher::operator()(const std::vector<int> vector) const {
    std::size_t seed = 0;
    // Iterates over vector values and combines their hashed values
    for (const int val : vector)
      boost::hash_combine(seed, boost::hash_value(val));
    return seed;
}

template std::size_t Hasher::operator()(const std::array<int, 1>) const;
template std::size_t Hasher::operator()(const std::array<int, 2>) const;
template std::size_t Hasher::operator()(const std::array<int, 3>) const;

} //namespace butil

#include <cstdlib>
#include <unordered_map>
#include <vector>

// Free functions that are helpful for writing tests

namespace btest {

//! Generates a random double between min and max
double RandomDouble(double min, double max) {
  double random_zero_one =
      static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  // multiply by total range and add to min
  return (random_zero_one * (max - min)) + min;
}

//! Generates a vector populated with n random doubles between min and max
std::vector<double> RandomVector(size_t n, double min, double max) {
  std::vector<double> return_vector;
  for (size_t i; i < n+1; ++i)
    return_vector.push_back(RandomDouble(min, max));
  return return_vector;
}

//! Generates a random unordered map of ints to vector<double>.
/*! Generates a random unordered map of ints to vector<double> map_size keys to
  vectors of length vector_size with values between min and max.
*/
std::unordered_map<int, std::vector<double>> RandomIntVectorMap(
    size_t map_size = 4, size_t vector_size = 3, double min = 0,
    double max = 100) {
  
  std::unordered_map<int, std::vector<double>> return_map;
  
  for (size_t i; i < map_size + 1; ++i) {
    int material_id = rand()%static_cast<int>(min - max + 1) + min;
    return_map.insert({material_id, RandomVector(vector_size, 0, 100)});
  }
  
  return return_map;
}

//! Overload to allow easier forward-declaration
std::unordered_map<int, std::vector<double>> RandomIntVectorMap() {
  return RandomIntVectorMap(4, 3, 0, 100);
}

} //namespace btest

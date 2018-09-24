#include <cstdlib>
#include <exception>
#include <unordered_map>
#include <vector>

#include <iostream>

#include <deal.II/lac/full_matrix.h>

// Free functions that are helpful for writing tests

namespace btest {

//! Generates a random double between min and max
double RandomDouble(double min, double max) {
  if (min >= max)
    throw std::runtime_error("Min must be less than max");
  double random_zero_one =
      static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  // multiply by total range and add to min
  return (random_zero_one * (max - min)) + min;
}

//! Generates a vector populated with n random doubles between min and max
std::vector<double> RandomVector(size_t n, double min, double max) {

  if (!n)
    throw std::runtime_error("Vector length must be > 0");
  
  std::vector<double> return_vector;
  for (size_t i = 0; i < n; ++i) {
    return_vector.push_back(RandomDouble(min, max));
  }
  return return_vector;
}

//! Generates a random unordered map of ints to vector<double>.
/*! Generates a random unordered map of ints to vector<double> map_size keys to
  vectors of length vector_size with values between min and max.
*/
std::unordered_map<int, std::vector<double>> RandomIntVectorMap(
    size_t map_size = 4, size_t vector_size = 3, double min = 0,
    double max = 100) {

  if (!map_size)
    throw std::runtime_error("IntVectorMap requires map size of at least 1");
  
  std::unordered_map<int, std::vector<double>> return_map;
  
  do {
    int material_id = rand()%(map_size*10);
    return_map.insert({material_id, RandomVector(vector_size, min, max)});
  } while (return_map.size() < map_size);
  
  return return_map;
}

//! Overload to allow easier forward-declaration
std::unordered_map<int, std::vector<double>> RandomIntVectorMap() {
  return RandomIntVectorMap(4, 3, 0, 100);
}

//! Generates a random dealii::FullMatrix<double>.
/*! Generates a random dealii::FullMatrix<double> of dimensions
  \f$\mathcal{R}^{m \times n}\f$ with random double values between min and max.
*/
  
dealii::FullMatrix<double> RandomMatrix(size_t m, size_t n, double min = 0,
                                        double max = 100) {
  dealii::FullMatrix<double> return_matrix(m, n);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      return_matrix.set(i, j, RandomDouble(min, max));
    }
  }                     
  return return_matrix;
}

//! Generates a random unordered map of ints to dealii::FullMatrix<double>
std::unordered_map<int, dealii::FullMatrix<double>> RandomIntMatrixMap(
    size_t map_size = 4, size_t m = 5, size_t n = 5, double min = 0,
    double max = 100) {
  
  std::unordered_map<int, dealii::FullMatrix<double>> return_map;
  
  for (size_t i; i < map_size + 1; ++i) {
    int material_id = rand()%static_cast<int>(min - max + 1) + min;
    return_map.insert({material_id, RandomMatrix(m, n, min, max)});
  }
  
  return return_map;
}

//! Overload for easy forward declaration
std::unordered_map<int, dealii::FullMatrix<double>> RandomIntMatrixMap() {
  return RandomIntMatrixMap(4, 5, 5, 0, 100);
}

} //namespace btest

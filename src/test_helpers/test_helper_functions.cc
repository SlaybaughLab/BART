#include "test_helper_functions.h"

// Free functions that are helpful for writing tests

namespace bart {

namespace test_helpers {

int RandomInt(int min, int max) {
  if (min >= max)
    throw std::runtime_error("Min must be less than max");
  double random_zero_one =
      static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
  return (random_zero_one * (max - min)) + min;
}

double RandomDouble(double min, double max) {
  if (min >= max)
    throw std::runtime_error("Min must be less than max");
  double random_zero_one =
      static_cast<double>(std::rand())/static_cast<double>(RAND_MAX);
  // multiply by total range and add to min
  return (random_zero_one * (max - min)) + min;
}

std::vector<double> RandomVector(std::size_t n, double min, double max) {

  if (n == 0u)
    throw std::runtime_error("Vector length must be > 0");

  std::vector<double> return_vector;
  for (std::size_t i = 0; i < n; ++i) {
    return_vector.push_back(RandomDouble(min, max));
  }
  return return_vector;
}

std::unordered_map<int, std::vector<double>>
RandomIntVectorMap(std::size_t map_size, std::size_t vector_size, double min,
                   double max) {
  if (map_size == 0u)
    throw std::runtime_error("IntVectorMap requires map size of at least 1");

  std::unordered_map<int, std::vector<double>> return_map;

  do {
    int material_id = std::rand()%(map_size*10);
    return_map.insert({material_id, RandomVector(vector_size, min, max)});
  } while (return_map.size() < map_size);

  return return_map;
}

dealii::FullMatrix<double> RandomMatrix(std::size_t m, std::size_t n, double min,
                                        double max) {
  dealii::FullMatrix<double> return_matrix(m, n);
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      return_matrix.set(i, j, RandomDouble(min, max));
    }
  }                     
  return return_matrix;
}

std::unordered_map<int, dealii::FullMatrix<double>>
RandomIntMatrixMap(std::size_t map_size, std::size_t m, std::size_t n,
                   double min, double max) {
  std::unordered_map<int, dealii::FullMatrix<double>> return_map;
  
  do {
    int material_id = std::rand()%(map_size*10);
    return_map.insert({material_id, RandomMatrix(m, n, min, max)});
  } while (return_map.size() < map_size);
  
  return return_map;
}

} // namespace test_helpers

} // namespace bart

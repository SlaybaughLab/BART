#ifndef BART_SRC_TEST_HELPERS_TEST_HELPER_FUNCTIONS_H_
#define BART_SRC_TEST_HELPERS_TEST_HELPER_FUNCTIONS_H_

#include <cstdlib>
#include <exception>
#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

namespace bart {

namespace test_helpers {

//! Generates a random double between min and max
double RandomDouble(double min, double max);

} // namespace test_helpers

} // namespace bart

namespace btest {

//! Generates a random double between min and max
[[deprecated]] double RandomDouble(double min, double max);

//! Generates a vector populated with n random doubles between min and max
std::vector<double> RandomVector(std::size_t n, double min, double max);

//! Generates a random unordered map of ints to vector<double>.
/*! Generates a random unordered map of ints to vector<double> map_size keys to
  vectors of length vector_size with values between min and max.
*/

std::unordered_map<int, std::vector<double>> RandomIntVectorMap(
    std::size_t map_size = 4, std::size_t vector_size = 3, double min = 0,
    double max = 100);

//! Generates a random dealii::FullMatrix<double>.
/*! Generates a random dealii::FullMatrix<double> of dimensions
  \f$\mathcal{R}^{m \times n}\f$ with random double values between min and max.
*/
dealii::FullMatrix<double> RandomMatrix(std::size_t m, std::size_t n,
                                        double min = 0, double max = 100);

//! Generates a random unordered map of ints to dealii::FullMatrix<double>

std::unordered_map<int, dealii::FullMatrix<double>>
RandomIntMatrixMap(std::size_t map_size = 4, std::size_t m = 5,
                   std::size_t n = 5, double min = 0, double max = 100);

} // namespace btest

#endif // BART_SRC_TEST_HELPERS_TEST_HELPER_FUNCTIONS_H_

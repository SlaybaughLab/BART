#ifndef BART_SRC_UTILITY_UTILITY_FUNCTIONS_HPP_
#define BART_SRC_UTILITY_UTILITY_FUNCTIONS_HPP_

#include <cmath>
#include <vector>

namespace bart::utility {

/*! \brief Conducts a precise summation by reducing truncation errors.
 * This precise summation is performed using the Neumaier variation of the Kahan summation.
 */
template <typename T>
auto PreciseSum(const T values) -> double {
  double sum = *values.begin();
  double correction = 0;
  for (auto cit = values.begin() + 1; cit < values.end(); ++cit) {
    double temporary_sum = sum + *cit;
    if (std::abs(sum) >= std::abs(*cit))
      correction += (sum - temporary_sum) + *cit;
    else
      correction += (*cit - temporary_sum) + sum;
    sum = temporary_sum;
  }
  return sum + correction;
}

} // namespace bart::utility

#endif // BART_SRC_UTILITY_UTILITY_FUNCTIONS_HPP_

#ifndef BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_
#define BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_

#include <cmath>
#include <vector>

namespace butil {

//! Conducts a precise summation by reducing truncation errors.
/*! Conduct a precise summation using the Neumaier variation of the Kahan
  summation.
*/
template <typename T>
double PreciseSum(const T values) {
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

} // namespace butil

#endif // BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_

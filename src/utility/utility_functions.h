#ifndef BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_
#define BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_

#include <vector>

namespace butil {

template <typename T>
double PreciseSum(const T values) {
  double sum = 0;
  double error = 0;
  for (const double& addend : values) {
    const double corrected_addend = addend - error;
    const double new_sum = sum + corrected_addend;
    // error is low order digits of corrected_addend lost when making new_sum
    error = (new_sum - sum) - corrected_addend;
    sum = new_sum;
  }
  return sum;
}

} // namespace butil

#endif // BART_SRC_UTILITY_UTILITY_FUNCTIONS_H_

#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_TESTS_K_EIGENVALUE_CALCULATOR_MOCK_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_TESTS_K_EIGENVALUE_CALCULATOR_MOCK_HPP_

#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::eigenvalue::k_eigenvalue {

class K_EigenvalueCalculatorMock : public K_EigenvalueCalculatorI {
 public:
  MOCK_METHOD(std::optional<double>, k_eigenvalue, (), (const, override));
  MOCK_METHOD(double, CalculateK_Eigenvalue, (system::System& system), (override));
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif //BART_SRC_EIGENVALUE_K_EIGENVALUE_TESTS_K_EIGENVALUE_CALCULATOR_MOCK_HPP_

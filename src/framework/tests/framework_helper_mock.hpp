#ifndef BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_HELPER_MOCK_HPP_
#define BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_HELPER_MOCK_HPP_

#include "framework/framework_helper_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::framework {

template <int dim>
class FrameworkHelperMock : public FrameworkHelperI<dim> {
 public:
  MOCK_METHOD(framework::FrameworkParameters, ToFrameworkParameters, (const problem::ParametersI& parameters), (override));
  MOCK_METHOD(std::unique_ptr<framework::FrameworkI>, BuildFramework, (builder::FrameworkBuilderI<dim>&,
      framework::FrameworkParameters&), (override));
  MOCK_METHOD(std::unique_ptr<framework::FrameworkI>, BuildFramework, (builder::FrameworkBuilderI<dim>&,
      framework::FrameworkParameters&, system::moments::SphericalHarmonicI*), (override));
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_HELPER_MOCK_HPP_

#ifndef BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_

#include "framework/builder/framework_builder_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::framework::builder {

template <int dim>
class FrameworkBuilderMock : public FrameworkBuilderI<dim> {
 public:
  using typename FrameworkBuilderI<dim>::FrameworkI;
  MOCK_METHOD(std::unique_ptr<FrameworkI>, BuildFramework, (const std::string&, const FrameworkParameters&), (override));
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_BUILDER_MOCK_HPP_

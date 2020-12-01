#ifndef BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_VALIDATOR_MOCK_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_VALIDATOR_MOCK_HPP_

#include "framework/builder/framework_validator_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::framework::builder {

class FrameworkValidatorMock : public FrameworkValidatorI {
 public:
  MOCK_METHOD(FrameworkValidatorI&, AddPart, (const FrameworkPart to_add), (override));
  MOCK_METHOD(void, Parse, (const framework::FrameworkParameters), (override));
  MOCK_METHOD(void, Parse, (const problem::ParametersI& to_parse), (override));
  MOCK_METHOD(void, ReportValidation, (), (override));

  MOCK_METHOD(bool, HasNeededParts, (), (override, const ));
  MOCK_METHOD(bool, HasUnneededParts, (), (override, const ));
  MOCK_METHOD(std::set<FrameworkPart>, NeededParts, (), (override, const ));
  MOCK_METHOD(std::set<FrameworkPart>, Parts, (), (override, const ));
  MOCK_METHOD(std::set<FrameworkPart>, UnneededParts, (), (override, const ));
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_TESTS_FRAMEWORK_VALIDATOR_MOCK_HPP_

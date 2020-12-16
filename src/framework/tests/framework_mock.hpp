#ifndef BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_MOCK_HPP_
#define BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_MOCK_HPP_

#include "framework/framework_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::framework {

class FrameworkMock : public FrameworkI {
  MOCK_METHOD(void, SolveSystem, (), (override));
  MOCK_METHOD(system::System*, system, (), (const, override));
  MOCK_METHOD(void, OutputResults, (std::ostream& output_stream), (override));
  MOCK_METHOD(void, OutputMasterFile, (std::ostream& output_stream, const std::vector<std::string>& filenames,
      const int process_id), (override));
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_TESTS_FRAMEWORK_MOCK_HPP_

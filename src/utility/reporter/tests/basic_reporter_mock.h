#ifndef BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_
#define BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_

#include "utility/reporter/basic_reporter_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace utility {

namespace reporter {

class BasicReporterMock : public BasicReporterI {
 public:
  MOCK_METHOD(void, Report, (const std::string &to_report), (override));
  MOCK_METHOD(void, Report, (const std::string &to_report, Color), (override));
  MOCK_METHOD(BasicReporterI&, Instream, (const std::string&));
  MOCK_METHOD(BasicReporterI&, Instream, (const Color));
  BasicReporterMock& operator<<(const std::string& to_report) override {
    Instream(to_report);
    return *this; }
  BasicReporterMock& operator<<(const Color color) override {
    Instream(color);
    return *this; }
};

} // namespace reporter

} // namespace utility

} // namespace bart

#endif //BART_SRC_UTILITY_REPORTER_TESTS_MPI_MOCK_H_

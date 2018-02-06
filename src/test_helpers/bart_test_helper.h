#ifndef BART_TEST_HELPER_H_
#define BART_TEST_HELPER_H_

#include <string>

namespace btest {

class BartTestHelper {
 public:
  BartTestHelper(bool report = false,
                 std::string gold_files_directory="/test_data");
  const std::string& GetReportDirectory() const { return report_directory_; };
 private:
  void MakeReportDirectory();
  const bool report_;
  const std::string gold_files_directory_;
  std::string report_directory_;
};

}

#endif

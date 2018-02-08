#ifndef BART_TEST_HELPER_H_
#define BART_TEST_HELPER_H_

#include <string>

namespace btest {

class BartTestHelper {
 public:
  BartTestHelper();
  BartTestHelper(bool report, std::string gold_files_directory);
  void ReInit(bool report, std::string gold_files_directory);
  const std::string& GetReportDirectory() const { return report_directory_; };
  bool GoldTest(std::string filename) const;
  void CleanupGold(std::string filename, bool result, bool actual_good) const;
  void MakeDiff(std::string filename, std::string diff) const;
 private:
  void MakeReportDirectory();
  bool report_;
  std::string gold_files_directory_;
  std::string report_directory_;
};
// This is a Global Bart Test Helper used for tests
BartTestHelper& GlobalBartTestHelper();

}

#endif

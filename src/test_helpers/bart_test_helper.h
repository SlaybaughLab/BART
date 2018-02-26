#ifndef BART_TEST_HELPER_H_
#define BART_TEST_HELPER_H_

#include <string>
#include <memory>
#include <fstream>
#include <sys/stat.h>
#include <ctime>
#include <iomanip>

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include "gtest/gtest.h"

#include "gold_stream_evaluator.h"

namespace btest {

class BartTestHelper {
 public:
  BartTestHelper();
  BartTestHelper(bool report, std::string gold_files_directory);
  
  void ReInit(bool report, std::string gold_files_directory);
  bool GoldTest(std::string filename) const;
  const std::string& GetReportDirectory() const { return report_directory_; };
  std::string GetFailMessage() const { return fail_message_; };
  void SetReport(bool report);
  void SetGoldFilesDirectory(std::string gold_files_directory);
  
  void OpenLog(std::string filename);
  void OpenMPILog(std::string filename);
  void CloseLog() { log_stream_.reset();};
  bool IsLogging() const { return log_stream_ != nullptr; };
  
 private:
  void CleanupGold(std::string filename, bool result, bool actual_good) const;
  void MakeDiff(std::string filename, std::string diff) const;
  void MakeReportDirectory();
  bool report_;
  std::string gold_files_directory_;
  std::string report_directory_;
  std::unique_ptr<std::ofstream> log_stream_;
  mutable std::string fail_message_ = "";
};
// This is a Global Bart Test Helper used for tests
BartTestHelper& GlobalBartTestHelper();

// Non-member helper functions
void GoldTestInit(std::string filename);
void GoldTestRun(std::string filename);
}

#endif

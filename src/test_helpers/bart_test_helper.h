#ifndef BART_TEST_HELPERS_BART_TEST_HELPER_H_
#define BART_TEST_HELPERS_BART_TEST_HELPER_H_

#include <ctime>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <sys/stat.h>

#include <deal.II/base/logstream.h>
#include <gtest/gtest.h>
#include "gold_stream_evaluator.h"

namespace btest {
//! This class provides a test framework for BART using gTest.
/*!

Functionality includes:

(1) Managing a filestream resource for the dealii log.

(2) Running a gold file test to compare two files.

 \author Joshua Rehak
 \date 2018/2
 */
class BARTTestHelper {
 public:
  //! Default constructor, default values `report = true`, `gold_files_directory = "test_data/"`
  BARTTestHelper();

  //! Constructor specifying if a report will be generated and location of gold files.
  /*!
    \param report a \ref bool indicating if a report will be generated. A value
    of `true` will immediately generate a report directory.
    \param gold_files_directory a std::string specifying the location of the
    gold files to be tested in the form `"dir/"`
  */
  BARTTestHelper(bool report, std::string gold_files_directory);

  //! Reinitializes the class
  void ReInit(bool report, std::string gold_files_directory);

  //! Runs a gold test, comparing `./filename` to `gold_files_directory/filename.gold`
  bool GoldTest(std::string filename) const;

  //! Returns the name of the report directory
  std::string GetReportDirectory() const;

  //! Returns a fail message if a gold test fails, or ""
  std::string GetFailMessage() const;

  //! Set the value of `report`
  void SetReport(bool report);

  //! Set the location of the gold files
  void SetGoldFilesDirectory(std::string gold_files_directory);

  //! Open attach an ofstream for ./filename to the deallii log
  void OpenLog(std::string filename);

  //! Detaches and closes the ofstream attached to the dealiilog
  void CloseLog();

  //! Returns the status of the dealii log
  bool IsLogging() const;

 private:
  //! Used to clean up files following a gold test
  void CleanupGold(std::string filename, bool result, bool actual_good) const;

  //! Generates the diff file between two files in unified format
  void MakeDiff(std::string filename, std::string diff) const;

  //! Creates the report directory
  void MakeReportDirectory();

  bool report_;

  std::string gold_files_directory_;

  std::string report_directory_;

  //! Stream that will be attached to the deallii log
  std::unique_ptr<std::ofstream> log_stream_;

  mutable std::string fail_message_ = "";
};

inline std::string BARTTestHelper::GetReportDirectory() const {
  return report_directory_;
}

inline std::string BARTTestHelper::GetFailMessage() const {
  return fail_message_;
}

inline void BARTTestHelper::CloseLog() {
  log_stream_.reset();
}

inline bool BARTTestHelper::IsLogging() const {
  return log_stream_ != nullptr;
}

/*! \relates BARTTestHelper
  This is the localstatic BARTTestHelper that is used by all gold file comparison tests
*/

BARTTestHelper& GlobalBARTTestHelper();

// Non-member helper functions
/*! \relates BARTTestHelper
  Use the global BARTTestHelper to initialize a gold test for a file,
  opens the dealii log outputting to `filename`
*/
void GoldTestInit(std::string filename);
/*! \relates BARTTestHelper
  Run the gold test for a file `filename`. This function closes the current
  stream to `filename`, and then compares it to `filename.gold` and then cleans up.
*/
void GoldTestRun(std::string filename);

//! This class provides a parallel gtest environment for BART.
/*!
 The main functionality is to provide MPI initilization and a tear down function.

 \author Weixiong Zheng
 \date 2018/8
 */
class BARTParallelEnvironment : public ::testing::Test {
 public:
  BARTParallelEnvironment();
  virtual ~BARTParallelEnvironment();

  virtual void TearDown();
 protected:
  void MPIInit();
  void MPIFinalize();
};
} // namespace btest

#endif // BART_TEST_HELPERS_BART_TEST_HELPER_H_

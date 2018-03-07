#ifndef BART_TEST_HELPERS_BART_TEST_HELPER_H_
#define BART_TEST_HELPERS_BART_TEST_HELPER_H_

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
//! This class provides a test framework for BART using gTest.
/*!
 
Functionality includes:

(1) Managing a filestream resource for the dealii log.

(2) Running a gold file test to compare two files.
 
 \author Joshua Rehak
 \date 2018/2
 */
class BartTestHelper {
 public:
  //! Default constructor, default values `report = true`, `gold_files_directory = "test_data/"`
  BartTestHelper();
  //! Constructor specifying if a report will be generated and location of gold files.
  /*!
    \param report a \ref bool indicating if a report will be generated. A value
    of `true` will immediately generate a report directory.
    \param gold_files_directory a std::string specifying the location of the
    gold files to be tested in the form `"dir/"`
  */
  BartTestHelper(bool report, std::string gold_files_directory);
  //! Reinitializes the class
  void ReInit(bool report, std::string gold_files_directory);
  //! Runs a gold test, comparing `./filename` to `gold_files_directory/filename.gold`
  bool GoldTest(std::string filename) const;
  //! Returns the name of the report directory
  const std::string& GetReportDirectory() const { return report_directory_; };
  //! Returns a fail message if a gold test fails, or ""
  std::string GetFailMessage() const { return fail_message_; };
  //! Set the value of `report`
  void SetReport(bool report);
  //! Set the location of the gold files
  void SetGoldFilesDirectory(std::string gold_files_directory);
  //! Open attach an ofstream for ./filename to the deallii log
  void OpenLog(std::string filename);
  //! Detaches and closes the ofstream attached to the dealiilog
  void CloseLog() { log_stream_.reset();};
  //! Returns the status of the dealii log
  bool IsLogging() const { return log_stream_ != nullptr; };
  
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
/*! \relates BartTestHelper
  This is the localstatic BartTestHelper that is used by all gold file comparison tests
*/

BartTestHelper& GlobalBartTestHelper();

// Non-member helper functions
/*! \relates BartTestHelper
  Use the global BartTestHelper to initialize a gold test for a file,
  opens the dealii log outputting to `filename`
*/
void GoldTestInit(std::string filename);
/*! \relates BartTestHelper
  Run the gold test for a file `filename`. This function closes the current
  stream to `filename`, and then compares it to `filename.gold` and then cleans up.
*/
void GoldTestRun(std::string filename);
} // namespace btest

#endif // BART_TEST_HELPERS_BART_TEST_HELPER_H_

#include "bart_test_helper.h"
#include <sys/stat.h>
#include <ctime>
#include <iomanip>
#include <fstream>

namespace btest {

BartTestHelper::BartTestHelper(bool report, std::string gold_files_directory)
    : report_(report),
      gold_files_directory_(gold_files_directory),
      report_directory_("") {
  if (report)
    MakeReportDirectory();
}

void BartTestHelper::MakeReportDirectory() {
  //Generate directory name
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::string directory_name;

  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M");

  std::string time = oss.str();

  // Name and create directory for failed tests
  report_directory_ = gold_files_directory_ + time + "_fail";

  //Make directory
  const int dir_err = mkdir(report_directory_.c_str(),
                              S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (-1 == dir_err)
    throw std::runtime_error(("Failed to create report directory " +
                              report_directory_).c_str());
}

}

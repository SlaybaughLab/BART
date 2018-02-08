#include "bart_test_helper.h"
#include "gold_stream_evaluator.h"
#include <sys/stat.h>
#include <memory>
#include <ctime>
#include <iomanip>
#include <fstream>

namespace btest {

BartTestHelper::BartTestHelper()
    : BartTestHelper(false, "test_data/") {}

BartTestHelper::BartTestHelper(bool report, std::string gold_files_directory)
    : report_(report),
      gold_files_directory_(gold_files_directory),
      report_directory_("") {
  if (report)
     MakeReportDirectory();
}

void BartTestHelper::ReInit(bool report, std::string gold_files_directory) {
  report_ = report;
  gold_files_directory_ = gold_files_directory;
  report_directory_ = "";
  if (report)
    MakeReportDirectory();
}

bool BartTestHelper::GoldTest(std::string filename) const {
  auto actual_file_stream = std::make_unique<std::ifstream>(filename);
  auto gold_file_stream =
      std::make_unique<std::ifstream>(gold_files_directory_ + filename + ".gold");
  
  GoldStreamEvaluator evaluator(std::move(gold_file_stream),
                                std::move(actual_file_stream));
  bool result = evaluator.RunGoldTest();
  
  // Make diff if required and possible
  if (!result && evaluator.ActualGood() && evaluator.GoldGood())
    MakeDiff(filename, evaluator.GetDiff());
    
  evaluator.CloseStreams();
  CleanupGold(filename, result, evaluator.ActualGood());
  
  // Messages about test failure
  if (!evaluator.GoldGood())
    printf("BartTestHelper: GoldTest Failed: Bad gold file\n");
  else if(!evaluator.ActualGood())
    printf("BartTestHelper: GoldTest Failed: No actual file\n");
  else if(!result)
    printf("BartTestHelper: GoldTest Failed, Actual/Gold not identical\n");

  return result;
}

void BartTestHelper::CleanupGold(std::string filename,
                                 bool result, bool actual_good) const {
  if (actual_good && (result || !report_)) {
    // Delete file
    const int remove_err = remove(filename.c_str());
    if (remove_err != 0)
      throw std::runtime_error(("Failed to delete actual test file: " +
                                filename).c_str());
  } else if (actual_good) {
    const int rename_err = rename(filename.c_str(),
                                  (report_directory_ + "/" + filename).c_str());
    if (rename_err != 0)
      throw std::runtime_error(("Failed to move actual test file: " +
                                filename + " to " +
                                report_directory_ + filename).c_str());
  }
}

void BartTestHelper::MakeDiff(std::string filename, std::string diff) const{
  std::string diff_filename = report_directory_ + "/" + filename + ".diff";
  std::ofstream diff_stream(diff_filename, std::ios_base::out);
  diff_stream << diff;
  diff_stream.close();
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

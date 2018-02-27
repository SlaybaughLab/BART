#include "bart_test_helper.h"

namespace btest {

BartTestHelper::BartTestHelper()
    : BartTestHelper(false, "test_data/") {}

BartTestHelper::BartTestHelper(bool report, std::string gold_files_directory) {
  ReInit(report, gold_files_directory);
}

void BartTestHelper::SetReport(bool report) {
  ReInit(report, gold_files_directory_);
}

void BartTestHelper::SetGoldFilesDirectory(std::string gold_files_directory) {
  ReInit(report_, gold_files_directory);
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

  // Make failure message
  if (!result) {
    if (!evaluator.ActualGood()) {
      fail_message_ = "Bad test file";
    } else if (!evaluator.GoldGood()) {
      fail_message_ = "Bad gold file";
    }
  }
  
  evaluator.CloseStreams();
  CleanupGold(filename, result, evaluator.ActualGood());
  
  return result;

}

void BartTestHelper::OpenLog(std::string filename) {
  if (log_stream_ == nullptr) {
    log_stream_ = std::make_unique<std::ofstream>(filename);
    dealii::deallog.attach(*log_stream_, false);
  } else
    throw std::runtime_error("BartTestHelper: Log is already open");
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
  //Generate directory if it isn't already made
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::string directory_name;

  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M_%S");

  std::string time = oss.str();

  // Name and create directory for failed tests
  report_directory_ = gold_files_directory_ + time + "_report";

  //Make directory
  struct stat sb;
  if (!(stat(report_directory_.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
    const int dir_err = mkdir(report_directory_.c_str(),
                              S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
      throw std::runtime_error(("Failed to create report directory " +
                                report_directory_).c_str());
  }
}

BartTestHelper& GlobalBartTestHelper() {
  static BartTestHelper global_bth;
  return global_bth;
}

void GoldTestInit(std::string filename) {
  GlobalBartTestHelper().OpenLog(filename);
}

void GoldTestRun(std::string filename) {
  GlobalBartTestHelper().CloseLog();
  ASSERT_TRUE(GlobalBartTestHelper().GoldTest(filename)) <<
      GlobalBartTestHelper().GetFailMessage();
}

}

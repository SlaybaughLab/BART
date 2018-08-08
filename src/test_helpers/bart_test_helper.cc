#include "bart_test_helper.h"

#include "gtest/gtest.h"

namespace btest {

BARTTestHelper::BARTTestHelper()
    : BARTTestHelper(false, "test_data/") {}

BARTTestHelper::BARTTestHelper(bool report, std::string gold_files_directory) {
  ReInit(report, gold_files_directory);
}

void BARTTestHelper::SetReport(bool report) {
  ReInit(report, gold_files_directory_);
}

void BARTTestHelper::SetGoldFilesDirectory(std::string gold_files_directory) {
  ReInit(report_, gold_files_directory);
}

void BARTTestHelper::ReInit(bool report, std::string gold_files_directory) {
  report_ = report;
  gold_files_directory_ = gold_files_directory;
  report_directory_ = "";
  if (report)
    MakeReportDirectory();
}

bool BARTTestHelper::GoldTest(std::string filename) const {
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
  //Close streams and clean up files
  evaluator.CloseStreams();
  CleanupGold(filename, result, evaluator.ActualGood());

  return result;

}

void BARTTestHelper::OpenLog(std::string filename) {
  // Make an ofstream targeting filename and attach to the dealii log, or return
  // an error if a log is already open
  if (log_stream_ == nullptr) {
    log_stream_ = std::make_unique<std::ofstream>(filename);
    dealii::deallog.attach(*log_stream_, false);
  } else
    throw std::runtime_error("BARTTestHelper: Log is already open");
}

void BARTTestHelper::CleanupGold(std::string filename,
                                 bool result, bool actual_good) const {
  // If the actual file is good, and the test passes or there is no report to
  // generate, delete the actual file
  if (actual_good && (result || !report_)) {
    const int remove_err = remove(filename.c_str());
    if (remove_err != 0)
      throw std::runtime_error(("Failed to delete actual test file: " +
                                filename).c_str());
  } else if (actual_good) {
    // Otherwise, if the actual file is good, move it to the report directory
    const int rename_err = rename(filename.c_str(),
                                  (report_directory_ + "/" + filename).c_str());
    if (rename_err != 0)
      throw std::runtime_error(("Failed to move actual test file: " +
                                filename + " to " +
                                report_directory_ + filename).c_str());
  }
}

void BARTTestHelper::MakeDiff(std::string filename, std::string diff) const{
  // Generate the diff file
  std::string diff_filename = report_directory_ + "/" + filename + ".diff";
  std::ofstream diff_stream(diff_filename, std::ios_base::out);
  diff_stream << diff;
  diff_stream.close();
}

void BARTTestHelper::MakeReportDirectory() {
  //Generate directory if it isn't already made
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::string directory_name;

  // Get the current time and date
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

BARTTestHelper& GlobalBARTTestHelper() {
  static BARTTestHelper global_bth;
  return global_bth;
}

void GoldTestInit(std::string filename) {
  GlobalBARTTestHelper().OpenLog(filename);
}

void GoldTestRun(std::string filename) {
  GlobalBARTTestHelper().CloseLog();
  ASSERT_TRUE(GlobalBARTTestHelper().GoldTest(filename)) <<
      GlobalBARTTestHelper().GetFailMessage();
}

BARTParallelEnvironment::BARTParallelEnvironment()
    : ::testing::Test() {}

BARTParallelEnvironment::~BARTParallelEnvironment() {}

void BARTParallelEnvironment::TearDown() {
  int err = MPI_Finalize();
  ASSERT_FALSE(err);
}

void BARTParallelEnvironment::MPIInit() {
  char** argv;
  int argc = 0;
  int err = MPI_Init (&argc, &argv);
  ASSERT_FALSE(err);
}

} // namespace btest

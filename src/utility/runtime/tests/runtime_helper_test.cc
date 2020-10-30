#include "utility/runtime/runtime_helper.h"

#include <sstream>
#include <vector>

#include <deal.II/base/mpi.h>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using ::testing::HasSubstr;

class UtilityRuntimeHelperTest : public ::testing::Test {
 public:
  const std::string version_{"7.6.9"}, filename{"testfile.input"};
  utility::runtime::RuntimeHelper test_helper{version_};

  void MakeArgv(std::string to_convert);

  std::vector<std::string> arguments_;
  std::vector<char*> argv_vector_;
  char** argv_;
  int argc_;
};

void UtilityRuntimeHelperTest::MakeArgv(std::string to_convert) {
  std::istringstream iss(to_convert);
  std::vector<std::string> arguments(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());
  arguments_ = arguments;
  for (auto& arg : arguments_) {
    argv_vector_.push_back(const_cast<char*>(arg.data()));
  }
  argv_vector_.push_back(nullptr);

  argv_ = const_cast<char**>(argv_vector_.data());
  argc_ = argv_vector_.size() - 1;
}

TEST_F(UtilityRuntimeHelperTest, Version) {
  EXPECT_EQ(test_helper.version(), version_);
}

TEST_F(UtilityRuntimeHelperTest, ProgramHeader) {
  auto header = test_helper.ProgramHeader();
  EXPECT_THAT(header, HasSubstr(version_));
}

TEST_F(UtilityRuntimeHelperTest, HelpMessage) {
  auto help_message = test_helper.HelpMessage();
  EXPECT_NE(help_message, "");
}

TEST_F(UtilityRuntimeHelperTest, PauseProgram) {
  MakeArgv("bart -p " + filename);
  EXPECT_FALSE(test_helper.do_pause());
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.do_pause());
  EXPECT_FALSE(test_helper.show_help());
  EXPECT_EQ(test_helper.filename(), filename);
}

TEST_F(UtilityRuntimeHelperTest, PauseProgramLong) {
  MakeArgv("bart --pause " + filename);
  EXPECT_FALSE(test_helper.do_pause());
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.do_pause());
  EXPECT_FALSE(test_helper.show_help());
  EXPECT_EQ(test_helper.filename(), filename);
}

TEST_F(UtilityRuntimeHelperTest, NoFileName) {
  MakeArgv("bart --pause");
  EXPECT_FALSE(test_helper.do_pause());
  EXPECT_FALSE(test_helper.show_help());
  EXPECT_ANY_THROW({
                     test_helper.ParseArguments(argc_, argv_);
                   });
}

TEST_F(UtilityRuntimeHelperTest, HelpLong) {
  MakeArgv("bart --help");
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.show_help());
}

TEST_F(UtilityRuntimeHelperTest, Help) {
  MakeArgv("bart -h");
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.show_help());
}

TEST_F(UtilityRuntimeHelperTest, UnknownOption) {
  MakeArgv("bart -k " + filename);
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.show_help());
}

TEST_F(UtilityRuntimeHelperTest, UnknownLongOption) {
  MakeArgv("bart --kill " + filename);
  EXPECT_FALSE(test_helper.show_help());
  test_helper.ParseArguments(argc_, argv_);
  EXPECT_TRUE(test_helper.show_help());
}




} // namespace

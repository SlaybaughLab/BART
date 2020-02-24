#include "test_helpers/gold_stream_evaluator.h"

#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class GoldStreamEvaluatorTest : public ::testing::Test {
 protected:
  virtual void SetUp();

  std::unique_ptr<std::istringstream> gold_iss;
  std::unique_ptr<std::istringstream> actual_iss;
};

void GoldStreamEvaluatorTest::SetUp() {
  gold_iss = std::make_unique<std::istringstream>();
  actual_iss = std::make_unique<std::istringstream>();
}

TEST_F(GoldStreamEvaluatorTest, BadGoldStream) {
  actual_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.GoldGood());
  ASSERT_TRUE(test_eval.ActualGood());
}

TEST_F(GoldStreamEvaluatorTest, BadActualStream) {
  gold_iss->setstate(std::ios_base::goodbit);
  actual_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.ActualGood());
  ASSERT_TRUE(test_eval.GoldGood());
}

TEST_F(GoldStreamEvaluatorTest, BadBothStream) {
  gold_iss->setstate(std::ios_base::badbit);
  actual_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.ActualGood());
  ASSERT_FALSE(test_eval.GoldGood());
}

TEST_F(GoldStreamEvaluatorTest, SameStream) {
  std::string input_text = "1\n2\n3\n4\n5";
  gold_iss->str(input_text);
  actual_iss->str(input_text);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_TRUE(test_eval.Compare());
}

TEST_F(GoldStreamEvaluatorTest, DiffStream) {
  std::string gold_text = "1\n2\n3\n4\n5";
  std::string actual_text = "1\n2\nX\n4\n5";
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.Compare());
}

TEST_F(GoldStreamEvaluatorTest, BadGoldStreamCompare) {
  actual_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_THROW(test_eval.Compare(), std::runtime_error);
}

TEST_F(GoldStreamEvaluatorTest, BadActualStreamCompare) {
  gold_iss->setstate(std::ios_base::goodbit);
  actual_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_THROW(test_eval.Compare(), std::runtime_error);
}

TEST_F(GoldStreamEvaluatorTest, CompareLongerGold) {
  std::string gold_text = "1\n2\n3\n4\n5";
  std::string actual_text = "1\n2\n3\n4";
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.Compare());
}

TEST_F(GoldStreamEvaluatorTest, CompareLongerActual) {
  std::string gold_text = "1\n2\n3\n4";
  std::string actual_text = "1\n2\n3\n4\n5";
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.Compare());
}

TEST_F(GoldStreamEvaluatorTest, DiffWorks) {
  std::string gold_text = "1\n2\nX\n4";
  std::string actual_text = "1\n2\n3\n4";
  std::string diff;
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);

  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  diff = test_eval.GetDiff();
  ASSERT_EQ(diff, "@@ -1,4 +1,4 @@\n 1\n 2\n-X\n+3\n 4\n");
}

TEST_F(GoldStreamEvaluatorTest, DiffWorksSame) {
  std::string gold_text = "1\n2\n3\n4";
  std::string actual_text = "1\n2\n3\n4";
  std::string diff;
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);

  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  diff = test_eval.GetDiff();
  ASSERT_EQ(diff, "");
}

TEST_F(GoldStreamEvaluatorTest, BadGoldStreamDiff) {
  actual_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_THROW(test_eval.GetDiff(), std::runtime_error);
}

TEST_F(GoldStreamEvaluatorTest, BadActualStreamDiff) {
  gold_iss->setstate(std::ios_base::goodbit);
  actual_iss->setstate(std::ios_base::badbit);
  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_THROW(test_eval.GetDiff(), std::runtime_error);
}

TEST_F(GoldStreamEvaluatorTest, RunGoldTestFail) {
  std::string gold_text = "1\n2\nX\n4";
  std::string actual_text = "1\n2\n3\n4";
  std::string diff;
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);

  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_FALSE(test_eval.RunGoldTest());
}

TEST_F(GoldStreamEvaluatorTest, RunGoldTestPass) {
  std::string gold_text = "1\n2\n3\n4";
  std::string actual_text = "1\n2\n3\n4";
  std::string diff;
  gold_iss->str(gold_text);
  actual_iss->str(actual_text);

  test_helpers::GoldStreamEvaluator test_eval(std::move(gold_iss),
                                       std::move(actual_iss));
  ASSERT_TRUE(test_eval.RunGoldTest());
}

} // namespace
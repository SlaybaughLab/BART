#include "gtest/gtest.h"

#include <memory>
#include <iostream>
#include <string>

#include "../stream_evaluator.h"

class StreamEvaluatorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    gold_iss = std::make_unique<std::istringstream>();
    temp_iss = std::make_unique<std::istringstream>();
  }

  std::unique_ptr<std::istringstream> gold_iss;
  std::unique_ptr<std::istringstream> temp_iss;

  btest::StreamEvaluator test_eval;
};

TEST_F(StreamEvaluatorTest, BadGoldStream)
{
  temp_iss->setstate(std::ios_base::goodbit);
  gold_iss->setstate(std::ios_base::badbit);
  test_eval.AdoptStreams(std::move(gold_iss), std::move(temp_iss));
  ASSERT_FALSE(test_eval.GoldGood());
  ASSERT_TRUE(test_eval.TempGood());
}

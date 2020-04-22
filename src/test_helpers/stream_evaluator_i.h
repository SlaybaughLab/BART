#ifndef BART_SRC_TEST_HELPERS_STREAM_EVALUATOR_I_H_
#define BART_SRC_TEST_HELPERS_STREAM_EVALUATOR_I_H_

#include <string>

namespace bart {

namespace test_helpers {

//! This class provides an interface for a class to evaluate two streams
/*!
  The two streams are the "gold" standard `gold_stream` and an "actual" stream
  `actual_stream`.
  \author Joshua Rehak
  \date 2018/2
*/
class StreamEvaluatorI {
 public:
  virtual ~StreamEvaluatorI() = default;
  
  //! Returns the result of the comparison
  virtual bool Compare() const = 0;
  
  //! Returns the difference between the two streams
  virtual std::string GetDiff() const = 0;
  
  //! Returns the result of a gold test on the two stream
  virtual bool RunGoldTest() const = 0;
  
  //! Returns the status of the gold stream
  virtual bool GoldGood() const = 0;
  
  //! Returns the status of the actual stream
  virtual bool ActualGood() const = 0;
};

} // namespace test_helpers

} // namespace bart

#endif // BART_SRC_TEST_HELPERS_STREAM_EVALUATOR_I_H_

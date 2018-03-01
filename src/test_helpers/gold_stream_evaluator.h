#ifndef GOLD_STREAM_EVALUATOR_H_
#define GOLD_STREAM_EVALUATOR_H_

#include <memory>
#include <iostream>

#include "stream_evaluator_i.h"

namespace btest {
//! This class will conduct a line by line comparison test on two streams.
/*!
 The two streams are the "gold" standard `gold_stream` and an "actual" stream
 `actual_stream`.
 
 \author Joshua Rehak
 \date 2018/2
 */
class GoldStreamEvaluator : public StreamEvaluatorI {
 public:
  //! Constructor, takes ownership of two streams for comparison
  GoldStreamEvaluator(std::unique_ptr<std::istream> gold_stream,
                    std::unique_ptr<std::istream> actual_stream);
  ~GoldStreamEvaluator() {};
  //! Returns a bool indicating if the streams are line-by-line identical
  bool Compare() const;
  //! Returns a diff string in unified format if the streams are different
  std::string GetDiff() const;
  /*! Returns the result of a "Gold" test, `true` if both streams are good
    and identical, `false` otherwise
  */
  bool RunGoldTest() const;
  //! Returns the status of the gold_stream
  const bool& GoldGood() const { return gold_good_; };
  //! Returns the status of the actual_stream
  const bool& ActualGood() const { return actual_good_; };
  //! Closes both stream
  void CloseStreams();
 private:
  //! Resets the streams to `bof`
  void ResetStreams() const;
  bool gold_good_ = false;
  bool actual_good_ = false;
  mutable std::unique_ptr<std::istream> gold_stream_;
  mutable std::unique_ptr<std::istream> actual_stream_;   
};

}

#endif

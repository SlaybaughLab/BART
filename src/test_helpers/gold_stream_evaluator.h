#ifndef GOLD_STREAM_EVALUATOR_H_
#define GOLD_STREAM_EVALUATOR_H_

#include <memory>
#include <iostream>

#include "stream_evaluator_i.h"

namespace btest {

class GoldStreamEvaluator : public StreamEvaluatorI {
 public:
  GoldStreamEvaluator(std::unique_ptr<std::istream> gold_stream,
                    std::unique_ptr<std::istream> actual_stream);
  ~GoldStreamEvaluator() {};
  bool Compare() const;
  std::string GetDiff() const;
  bool RunGoldTest() const;
  const bool& GoldGood() const { return gold_good_; };
  const bool& ActualGood() const { return actual_good_; };
  void CloseStreams();
 private:
  void ResetStreams() const;
  bool gold_good_ = false;
  bool actual_good_ = false;
  mutable std::unique_ptr<std::istream> gold_stream_;
  mutable std::unique_ptr<std::istream> actual_stream_;
    
};

}

#endif

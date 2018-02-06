#ifndef GOLD_STREAM_EVALUATOR_H_
#define GOLD_STREAM_EVALUATOR_H_

#include <memory>
#include <iostream>

#include "stream_evaluator_i.h"

namespace btest {

class GoldStreamEvaluator : public StreamEvaluatorI {
 public:
  GoldStreamEvaluator() {};
  ~GoldStreamEvaluator() {};
  void AdoptStreams(std::unique_ptr<std::istream> gold_stream,
                    std::unique_ptr<std::istream> actual_stream);
  bool Compare();
  std::string GetDiff();
  bool RunGoldTest();
  bool GoldGood() { return gold_good_; };
  bool ActualGood() { return actual_good_; };
   
 private:
  void ResetStreams();
  bool gold_good_ = false;
  bool actual_good_ = false;
  std::unique_ptr<std::istream> gold_stream_;
  std::unique_ptr<std::istream> actual_stream_;
    
};

}

#endif

#ifndef STREAM_EVALUATOR_H_
#define STREAM_EVALUATOR_H_

#include <memory>
#include <iostream>

#include "stream_evaluator_i.h"

namespace btest {

class StreamEvaluator : public StreamEvaluatorI {
 public:
  StreamEvaluator() {};
  ~StreamEvaluator() {};
  void AdoptStreams(std::unique_ptr<std::istream> gold_stream,
                    std::unique_ptr<std::istream> temp_stream);
  bool Compare();
  std::string GetDiff();
  bool GoldGood() { return gold_good_; };
  bool TempGood() { return temp_good_; };
   
 private:
  void ResetStreams();
  bool gold_good_ = false;
  bool temp_good_ = false;
  std::unique_ptr<std::istream> gold_stream_;
  std::unique_ptr<std::istream> temp_stream_;
    
};

}

#endif

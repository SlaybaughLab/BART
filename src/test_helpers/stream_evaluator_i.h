#ifndef STREAM_EVALUATOR_I_H_
#define STREAM_EVALUATOR_I_H_

#include <string>

namespace btest {

class StreamEvaluatorI {
 public:
  virtual ~StreamEvaluatorI() {};
  virtual bool Compare() = 0;
  virtual std::string GetDiff() = 0;
  virtual bool GoldGood() = 0;
  virtual bool ActualGood() = 0;
};


}

#endif

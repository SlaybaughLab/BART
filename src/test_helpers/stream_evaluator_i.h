#ifndef STREAM_EVALUATOR_I_H_
#define STREAM_EVALUATOR_I_H_

#include <string>

namespace btest {

class StreamEvaluatorI {
 public:
  virtual ~StreamEvaluatorI() {};
  virtual bool Compare() const = 0;
  virtual std::string GetDiff() const = 0;
  virtual bool RunGoldTest() const = 0;
  virtual bool GoldGood() const = 0;
  virtual bool ActualGood() const = 0;
};


}

#endif

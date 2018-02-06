#include "gmock/gmock.h"
#include "../stream_evaluator_i.h"

namespace btest {

class MockStreamEvaluator : public StreamEvaluatorI {
 public:
  MOCK_CONST_METHOD0(Compare, bool());
  MOCK_CONST_METHOD0(GetDiff, std::string());
  MOCK_CONST_METHOD0(RunGoldTest, bool());
  MOCK_CONST_METHOD0(GoldGood, const bool&());
  MOCK_CONST_METHOD0(ActualGood, const bool&());
};

}

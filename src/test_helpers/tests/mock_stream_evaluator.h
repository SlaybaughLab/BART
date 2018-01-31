#include "gmock/gmock.h"
#include "../stream_evaluator_i.h"

namespace btest {

class MockStreamEvaluator : public StreamEvaluatorI {
 public:
  MOCK_METHOD0(Compare, bool());
  MOCK_METHOD0(GetDiff, std::string());
  MOCK_METHOD0(GoldGood, bool());
  MOCK_METHOD0(TempGood, bool());
};

}

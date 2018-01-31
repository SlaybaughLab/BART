#include "stream_evaluator.h"

namespace btest {

void StreamEvaluator::AdoptStreams(std::unique_ptr<std::istream> gold_stream,
                                   std::unique_ptr<std::istream> temp_stream) {
  gold_stream_ = std::move(gold_stream);
  temp_stream_ = std::move(temp_stream);
  gold_good_ = gold_stream_->good();
  temp_good_ = temp_stream_->good();
}


}

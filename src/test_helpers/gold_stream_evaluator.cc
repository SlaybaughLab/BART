#include "gold_stream_evaluator.h"
#include "dtl/dtl.hpp"

#include <sstream>
#include <vector>

namespace btest {

void GoldStreamEvaluator::AdoptStreams(std::unique_ptr<std::istream> gold_stream,
                                   std::unique_ptr<std::istream> actual_stream) {
  gold_stream_ = std::move(gold_stream);
  actual_stream_ = std::move(actual_stream);
  gold_good_ = gold_stream_->good();
  actual_good_ = actual_stream_->good();
}

bool GoldStreamEvaluator::Compare() {
  if (!gold_good_ || !actual_good_)
    throw std::runtime_error("Cannot compare bad streams");

  ResetStreams();

  std::string gold_line;
  std::string actual_line;
  bool same = true;

  while (getline(*gold_stream_, gold_line))
  {
    //Check if files are identical
    getline(*actual_stream_, actual_line);
    if (gold_line != actual_line)
    {
      same = false;
    }
  }

  if (!gold_stream_->eof() || !actual_stream_->eof())
  {
    //One file longer than the other
    same = false;
    
  }
  
  return same;
}

std::string GoldStreamEvaluator::GetDiff() {
  if (!gold_good_ || !actual_good_)
    throw std::runtime_error("Cannot diff bad streams");
  
  //Generate diff and return output stream with it
  std::ostringstream diff_stream;
  std::string line_buffer;
  std::vector<std::string> gold_lines, actual_lines;

  ResetStreams();

  while(getline(*gold_stream_, line_buffer))
    gold_lines.push_back(line_buffer);
  while(getline(*actual_stream_, line_buffer))
    actual_lines.push_back(line_buffer);

  dtl::Diff<std::string> diff(gold_lines, actual_lines);
  diff.onHuge();
  diff.compose();
  diff.composeUnifiedHunks();
  diff.printUnifiedFormat(diff_stream);
  
  return diff_stream.str();
}

bool GoldStreamEvaluator::RunGoldTest() {
  if (!gold_good_ || !actual_good_)
    return false;
  return Compare();    
}

void GoldStreamEvaluator::ResetStreams() {
  gold_stream_->seekg(0, std::ios::beg);
  actual_stream_->seekg(0, std::ios::beg);
}

}

#include "stream_evaluator.h"
#include "dtl/dtl.hpp"

#include <sstream>
#include <vector>

namespace btest {

void StreamEvaluator::AdoptStreams(std::unique_ptr<std::istream> gold_stream,
                                   std::unique_ptr<std::istream> temp_stream) {
  gold_stream_ = std::move(gold_stream);
  temp_stream_ = std::move(temp_stream);
  gold_good_ = gold_stream_->good();
  temp_good_ = temp_stream_->good();
}

bool StreamEvaluator::Compare() {
  if (!gold_good_ || !temp_good_)
    throw std::runtime_error("Cannot compare bad streams");

  ResetStreams();

  std::string gold_line;
  std::string temp_line;
  bool same = true;

  while (getline(*gold_stream_, gold_line))
  {
    //Check if files are identical
    getline(*temp_stream_, temp_line);
    if (gold_line != temp_line)
    {
      same = false;
    }
  }

  if (!gold_stream_->eof() || !temp_stream_->eof())
  {
    //One file longer than the other
    same = false;
    
  }
  
  return same;
}

std::string StreamEvaluator::GetDiff() {
  if (!gold_good_ || !temp_good_)
    throw std::runtime_error("Cannot diff bad streams");
  
  //Generate diff and return output stream with it
  std::ostringstream diff_stream;
  std::string line_buffer;
  std::vector<std::string> gold_lines, temp_lines;

  ResetStreams();

  while(getline(*gold_stream_, line_buffer))
    gold_lines.push_back(line_buffer);
  while(getline(*temp_stream_, line_buffer))
    temp_lines.push_back(line_buffer);

  dtl::Diff<std::string> diff(gold_lines, temp_lines);
  diff.onHuge();
  diff.compose();
  diff.composeUnifiedHunks();
  diff.printUnifiedFormat(diff_stream);
  
  return diff_stream.str();
}

void StreamEvaluator::ResetStreams() {
  gold_stream_->seekg(0, std::ios::beg);
  temp_stream_->seekg(0, std::ios::beg);
}

}

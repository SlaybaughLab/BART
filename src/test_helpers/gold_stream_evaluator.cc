#include "gold_stream_evaluator.h"
#include "dtl/dtl.hpp"

#include <sstream>
#include <vector>

namespace bart {

namespace test_helpers {

GoldStreamEvaluator::GoldStreamEvaluator(
    std::unique_ptr<std::istream> gold_stream,
    std::unique_ptr<std::istream> actual_stream) {
  // Take ownership of the streams and verify they are good
  gold_stream_ = std::move(gold_stream);
  actual_stream_ = std::move(actual_stream);
  gold_good_ = gold_stream_->good();
  actual_good_ = actual_stream_->good();
}

bool GoldStreamEvaluator::Compare() const {
  if (!gold_good_ || !actual_good_)
    throw std::runtime_error("Cannot compare bad streams");

  ResetStreams();

  std::string gold_line;
  std::string actual_line;
  bool same = true;
  unsigned int c1=0, c2=0;

  // Count the lines of each stream
  while(!gold_stream_->eof())  {
    getline(*gold_stream_, gold_line);
    ++c1;
  }
  while(!actual_stream_->eof())  {
    getline(*actual_stream_, actual_line);
    ++c2;
  }
  //Return to the beginning of the stream
  ResetStreams();

  // Do a line-by-line comparison of the streams if they have the same number
  // of lines
  if (c1 != c2) {
    same = false;
  } else {
    while (getline(*gold_stream_, gold_line)) {
      //Check if files are identical
      getline(*actual_stream_, actual_line);
      if (gold_line != actual_line)
      {
        same = false;
      }
    }
  }

  return same;
}

std::string GoldStreamEvaluator::GetDiff() const {
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

bool GoldStreamEvaluator::RunGoldTest() const {
  if (!gold_good_ || !actual_good_)
    return false;
  return Compare();    
}


void GoldStreamEvaluator::CloseStreams() {
  // Release the streams
  gold_stream_.reset();
  actual_stream_.reset();
}

void GoldStreamEvaluator::ResetStreams() const {
  // Clear any eof flags from the streams and return to the start of the stream
  gold_stream_->clear();
  actual_stream_->clear();
  gold_stream_->seekg(0, std::ios::beg);
  actual_stream_->seekg(0, std::ios::beg);
}

} // namespace test_helpers

} // namespace bart

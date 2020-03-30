#ifndef BART_SRC_UTILITY_REPORTER_MPI_H_
#define BART_SRC_UTILITY_REPORTER_MPI_H_

#include <memory>

#include <deal.II/base/conditional_ostream.h>
#include <ostream>

#include "utility/reporter/basic_reporter_i.h"
#include "utility/uncopyable.h"

namespace bart {

namespace utility {

namespace reporter {

/*! \brief Noisy reporter, will report each time Report is called for
 * convergence, or a passed string.
 */

class Mpi : public BasicReporterI, private Uncopyable {
 public:
  Mpi(std::unique_ptr<dealii::ConditionalOStream> pout_ptr)
      : pout_ptr_(std::move(pout_ptr)) {};
  ~Mpi() = default;

  void Report(const std::string &to_report) override {
    *pout_ptr_ << to_report;
  }

  virtual void Report(const std::string &to_report, Color color) override {
    *pout_ptr_ << color_string_.at(color) + to_report + color_string_.at(Color::Reset);
  }

 private:
  std::unique_ptr<dealii::ConditionalOStream> pout_ptr_;
  std::unordered_map<Color, std::string> color_string_{
      {Color::Reset, "\033[0m"},
      {Color::Red,   "\033[31m"},
      {Color::Green, "\033[32m"},
      {Color::Blue,  "\033[34m"},
  };
};
} // namespace reporter

} // namespace utility

} // namespace bart

#endif // BART_SRC_UTILITY_REPORTER_MPI_H_
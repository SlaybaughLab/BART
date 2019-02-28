#ifndef BART_SRC_CONVERGERNCE_FINAL_I_H_
#define BART_SRC_CONVERGERNCE_FINAL_I_H_

#include <optional>

namespace bart {

namespace convergence {

/*! \brief Interface for classes that check for final convergence of a
 * framework.
 */

class FinalI {
  /*! Contains the status of a convergence check */
  struct Status {
    int iteration_number = 0;
    int max_iterations = 0;
    bool is_complete = false;
    std::optional<int> failed_index = std::nullopt;
    std::optional<double> delta = std::nullopt;
  };
  virtual ~FinalI() = default;
  virtual Status CheckFinalConvergence();
  virtual Status convergence_status() const;
  virtual bool   convergence_is_complete() const;
  virtual int    max_iterations() const;
  virtual void   SetMaxIterations(int to_set);
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGERNCE_FINAL_I_H_

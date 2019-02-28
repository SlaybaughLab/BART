#ifndef BART_SRC_CONVERGERNCE_FINAL_I_H_
#define BART_SRC_CONVERGERNCE_FINAL_I_H_

#include <optional>

#include "convergence/status.h"

namespace bart {

namespace convergence {

/*! \brief Interface for classes that check for final convergence of a
 * framework.
 */

class FinalI {
 public:
  virtual ~FinalI() = default;
  virtual Status CheckFinalConvergence() = 0;
  virtual Status convergence_status() const = 0;
  virtual bool   convergence_is_complete() const = 0;
  virtual int    max_iterations() const = 0;
  virtual void   SetMaxIterations(int to_set) = 0;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGERNCE_FINAL_I_H_

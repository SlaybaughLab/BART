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
  virtual Status CheckFinalConvergence();
  virtual Status convergence_status() const;
  virtual bool   convergence_is_complete() const;
  virtual int    max_iterations() const;
  virtual void   SetMaxIterations(int to_set);
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGERNCE_FINAL_I_H_

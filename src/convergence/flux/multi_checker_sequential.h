#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_SEQUENTIAL_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_SEQUENTIAL_H_

#include <memory>

#include "convergence/flux/single_checker_i.h"
#include "convergence/flux/multi_checker.h"
#include "data/vector_parameters.h"


namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks each flux group sequentially for convergence using a
 * provided flux checker. Stops after the first non-converged flux is found */

class MultiCheckerSequential : public MultiChecker {
 public:
  explicit MultiCheckerSequential(std::unique_ptr<SingleCheckerI> checker)
      : MultiChecker(std::move(checker)) {};
  ~MultiCheckerSequential() = default; 
  bool CheckIfConverged(data::MultiFluxPtrs &current_iteration,
                        data::MultiFluxPtrs &previous_iteration) override;
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_SEQUENTIAL_H_

#include "convergence/final.h"
#include "final.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace convergence {

void Final::SetMaxIterations(IterationNumber to_set) {
    AssertThrow(to_set > 0,
                dealii::ExcMessage("Max iterations must be > 0"));
    max_iterations_ = to_set;
}
void Final::SetIteration(IterationNumber to_set) {
    AssertThrow(to_set > 0,
                dealii::ExcMessage("Iteration must be > 0"));
    iteration_ = to_set;
}

} // namespace convergence

} // namespace bart

#include "convergence/final.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace convergence {

Final& Final::SetMaxIterations(IterationNumber to_set) {
    AssertThrow(to_set > 0,
                dealii::ExcMessage("Max iterations must be > 0"));
    max_iterations_ = to_set;
    return *this;
}
Final& Final::SetIteration(IterationNumber to_set) {
    AssertThrow(to_set >= 0,
                dealii::ExcMessage("Iteration must be >= 0"));
    iteration_ = to_set;
    return *this;
}

} // namespace convergence

} // namespace bart

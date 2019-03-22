#include "convergence/final.h"

#include "convergence/moments/single_moment_checker_i.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace convergence {

Final& Final::SetMaxIterations(
    IterationNumber to_set) {

    AssertThrow(to_set > 0,
                dealii::ExcMessage("Max iterations must be > 0"));
    convergence_status_.max_iterations = to_set;
    return *this;
}

Final& Final::SetIteration(
    IterationNumber to_set) {

    AssertThrow(to_set >= 0,
                dealii::ExcMessage("Iteration must be >= 0"));
    convergence_status_.iteration_number = to_set;
    return *this;
}


} // namespace convergence

} // namespace bart

#include "convergence/final.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace convergence {

void Final::SetMaxIterations(int to_set) {
    AssertThrow(to_set > 0,
                dealii::ExcMessage("Max iterations must be > 0"));
    max_iterations_ = to_set;
}               

} // namespace convergence

} // namespace bart

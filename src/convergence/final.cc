#include "convergence/final.h"

#include "system/moments/spherical_harmonic_types.h"

#include <deal.II/base/exceptions.h>

namespace bart {

namespace convergence {

template <typename CompareType>
Final<CompareType>& Final<CompareType>::SetMaxIterations(
    IterationNumber to_set) {

    AssertThrow(to_set > 0,
                dealii::ExcMessage("Max iterations must be > 0"));
    convergence_status_.max_iterations = to_set;
    return *this;
}

template <typename CompareType>
Final<CompareType>& Final<CompareType>::SetIteration(
    IterationNumber to_set) {

    AssertThrow(to_set >= 0,
                dealii::ExcMessage("Iteration must be >= 0"));
    convergence_status_.iteration_number = to_set;
    return *this;
}

template class Final<system::moments::MomentVector>;
template class Final<system::moments::MomentsMap>;
template class Final<double>;


} // namespace convergence

} // namespace bart

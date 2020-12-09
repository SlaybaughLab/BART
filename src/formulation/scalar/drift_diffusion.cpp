#include "formulation/scalar/drift_diffusion.hpp"

namespace bart::formulation::scalar {

template class DriftDiffusion<1>;
template class DriftDiffusion<2>;
template class DriftDiffusion<3>;

} // namespace bart::formulation::scalar

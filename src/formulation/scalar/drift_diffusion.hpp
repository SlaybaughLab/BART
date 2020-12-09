#ifndef BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_
#define BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_

#include "formulation/scalar/drift_diffusion_i.hpp"

namespace bart::formulation::scalar {

template <int dim>
class DriftDiffusion : public DriftDiffusionI<dim> {
 public:
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_HPP_

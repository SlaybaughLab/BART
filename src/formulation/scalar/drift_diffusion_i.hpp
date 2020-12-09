#ifndef BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_
#define BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_

namespace bart::formulation::scalar {

template <int dim>
class DriftDiffusionI {
 public:
  virtual ~DriftDiffusionI() = default;
};

} // namespace bart::formulation::scalar

#endif //BART_SRC_FORMULATION_SCALAR_DRIFT_DIFFUSION_I_HPP_

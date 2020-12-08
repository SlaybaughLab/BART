#ifndef BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_I_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_I_HPP_

#include <map>

#include <deal.II/lac/vector.h>

#include "quadrature/quadrature_types.h"

namespace bart::quadrature::calculators {

class DriftDiffusionIntegratedFluxI {
 public:
  virtual ~DriftDiffusionIntegratedFluxI() = default;
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<dealii::Vector<double>>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, VectorPtr>;

  virtual auto Integrate(const VectorMap&) const -> std::vector<double> = 0;
};

} // namespace bart::quadrature::calculators

#endif //BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_I_HPP_

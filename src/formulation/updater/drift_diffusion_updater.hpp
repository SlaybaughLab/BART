#ifndef BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_
#define BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_

#include "diffusion_updater.hpp"
#include "formulation/scalar/drift_diffusion_i.hpp"
#include "formulation/stamper_i.h"
#include "quadrature/calculators/drift_diffusion_integrated_flux_i.hpp"

namespace bart::formulation::updater {

template <int dim>
class DriftDiffusionUpdater : public DiffusionUpdater<dim> {
 public:
  using DiffusionFormulation = typename DiffusionUpdater<dim>::DiffusionFormulationType;
  using Stamper = typename DiffusionUpdater<dim>::StamperType;

  DriftDiffusionUpdater(std::unique_ptr<DiffusionFormulation>,
                        std::shared_ptr<Stamper>,
                        std::unordered_set<problem::Boundary> reflective_boundaries = {});
  virtual ~DriftDiffusionUpdater() = default;
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_

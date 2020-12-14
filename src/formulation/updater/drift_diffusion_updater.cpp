#include "drift_diffusion_updater.hpp"

namespace bart::formulation::updater {

template<int dim>
DriftDiffusionUpdater<dim>::DriftDiffusionUpdater(std::unique_ptr<DiffusionFormulation> diffusion_formulation_ptr,
                                                  std::shared_ptr<Stamper> stamper_ptr,
                                                  std::unordered_set<problem::Boundary> reflective_boundaries)
    : DiffusionUpdater<dim>(std::move(diffusion_formulation_ptr), stamper_ptr, reflective_boundaries) {}

} // namespace bart::formulation::updater

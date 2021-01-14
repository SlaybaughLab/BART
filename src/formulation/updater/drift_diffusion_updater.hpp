#ifndef BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_
#define BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_

#include "diffusion_updater.hpp"
#include "formulation/scalar/drift_diffusion_i.hpp"
#include "formulation/stamper_i.h"
#include "quadrature/calculators/angular_flux_integrator_i.hpp"
#include "system/solution/solution_types.h"
#include "system/moments/spherical_harmonic_i.h"
#include "utility/has_dependencies.h"

namespace bart::formulation::updater {

template <int, typename ...T> class DriftDiffusionUpdaterFactory;

template <int dim>
class DriftDiffusionUpdater : public DiffusionUpdater<dim>, public utility::HasDependencies {
 public:
  using AngularFluxStorageMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using DiffusionFormulation = typename DiffusionUpdater<dim>::DiffusionFormulationType;
  using DriftDiffusionFormulation = formulation::scalar::DriftDiffusionI<dim>;
  using IntegratedFluxCalculator = quadrature::calculators::AngularFluxIntegratorI;
  using HighOrderMoments = system::moments::SphericalHarmonicI;
  using Stamper = typename DiffusionUpdater<dim>::StamperType;

  using Factory = DriftDiffusionUpdaterFactory<dim, std::unique_ptr<DiffusionFormulation>, std::unique_ptr<DriftDiffusionFormulation>, std::shared_ptr<Stamper>, std::shared_ptr<IntegratedFluxCalculator>,std::shared_ptr<HighOrderMoments>, AngularFluxStorageMap&, std::unordered_set<problem::Boundary>>;

  DriftDiffusionUpdater(std::unique_ptr<DiffusionFormulation>,
                        std::unique_ptr<DriftDiffusionFormulation>,
                        std::shared_ptr<Stamper>,
                        std::shared_ptr<IntegratedFluxCalculator>,
                        std::shared_ptr<HighOrderMoments>,
                        AngularFluxStorageMap&,
                        std::unordered_set<problem::Boundary> reflective_boundaries = {});
  virtual ~DriftDiffusionUpdater() = default;


  auto angular_flux_storage_map() const -> AngularFluxStorageMap {
    return angular_flux_storage_map_; }
  auto drift_diffusion_formulation_ptr() const -> DriftDiffusionFormulation* {
    return drift_diffusion_formulation_ptr_.get(); }
  auto high_order_moments() const -> HighOrderMoments* {
    return high_order_moments_.get(); }
  auto integrated_flux_calculator_ptr() const -> IntegratedFluxCalculator* {
    return integrated_flux_calculator_ptr_.get(); }
 protected:
  auto SetUpFixedFunctions(system::System&, system::EnergyGroup, quadrature::QuadraturePointIndex) -> void override;
  AngularFluxStorageMap angular_flux_storage_map_{};
  std::shared_ptr<HighOrderMoments> high_order_moments_;
  std::unique_ptr<DriftDiffusionFormulation> drift_diffusion_formulation_ptr_{ nullptr };
  std::shared_ptr<IntegratedFluxCalculator> integrated_flux_calculator_ptr_{ nullptr };
 private:
  static bool is_registered_;
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_DRIFT_DIFFUSION_UPDATER_HPP_

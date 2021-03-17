#ifndef BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_

#include <memory>
#include <unordered_set>

#include "formulation/scalar/diffusion_i.hpp"
#include "formulation/stamper_i.hpp"
#include "formulation/updater/fixed_updater.hpp"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/fixed_source_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "quadrature/quadrature_set_i.hpp"
#include "problem/parameter_types.hpp"
#include "utility/has_description.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class DiffusionUpdater
    : public FixedUpdater<dim>, public ScatteringSourceUpdaterI,
      public FissionSourceUpdaterI, public FixedSourceUpdaterI,
      public utility::HasDescription {
 public:
  using typename FixedUpdater<dim>::CellPtr;
  using typename FixedUpdater<dim>::MatrixFunction;
  using typename FixedUpdater<dim>::VectorFunction;
  using typename FixedUpdater<dim>::MatrixBoundaryFunction;
  using typename FixedUpdater<dim>::VectorBoundaryFunction;

  using DiffusionFormulationType = formulation::scalar::DiffusionI<dim>;
  using StamperType = formulation::StamperI<dim>;
  DiffusionUpdater(std::unique_ptr<DiffusionFormulationType>,
                   std::shared_ptr<StamperType>,
                   std::unordered_set<problem::Boundary> reflective_boundaries = {});
  virtual ~DiffusionUpdater() = default;

  void UpdateScatteringSource(
      system::System &,
      system::EnergyGroup,
      quadrature::QuadraturePointIndex) override;

  void UpdateFissionSource(
      system::System &,
      system::EnergyGroup,
      quadrature::QuadraturePointIndex) override;

  void UpdateFixedSource(
      system::System &to_update,
      system::EnergyGroup group,
      quadrature::QuadraturePointIndex index) override;

  std::unordered_set<problem::Boundary>& reflective_boundaries() {
    return reflective_boundaries_; }
  DiffusionFormulationType* formulation_ptr() const {
    return formulation_ptr_.get(); }
  StamperType* stamper_ptr() const { return stamper_ptr_.get(); }
 protected:
  auto SetUpFixedFunctions(system::System&, system::EnergyGroup, quadrature::QuadraturePointIndex) -> void override;
  std::unique_ptr<DiffusionFormulationType> formulation_ptr_;
  std::shared_ptr<StamperType> stamper_ptr_;
  std::unordered_set<problem::Boundary> reflective_boundaries_;
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_DIFFUSION_UPDATER_H_

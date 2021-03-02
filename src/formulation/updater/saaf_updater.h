#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_

#include <memory>
#include <unordered_set>

#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/stamper_i.hpp"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/boundary_conditions_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "quadrature/quadrature_set_i.h"
#include "problem/parameter_types.hpp"
#include "system/solution/solution_types.h"
#include "utility/has_description.h"

namespace bart {

namespace formulation {

namespace updater {

template <int dim>
class SAAFUpdater :
    public FixedUpdaterI,
    public BoundaryConditionsUpdaterI,
    public ScatteringSourceUpdaterI,
    public FissionSourceUpdaterI,
    public utility::HasDescription {
 public:
  using Boundary = problem::Boundary;
  using EnergyGroupToAngularSolutionPtrMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using StamperType = formulation::StamperI<dim>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;

  SAAFUpdater(std::unique_ptr<SAAFFormulationType>,
              std::unique_ptr<StamperType>,
              const std::shared_ptr<QuadratureSetType>&);
  SAAFUpdater(std::unique_ptr<SAAFFormulationType>,
              std::unique_ptr<StamperType>,
              const std::shared_ptr<QuadratureSetType>&,
              const EnergyGroupToAngularSolutionPtrMap&,
              const std::unordered_set<Boundary>);

  void UpdateBoundaryConditions(system::System &to_update,
                                system::EnergyGroup group,
                                quadrature::QuadraturePointIndex index) override;
  void UpdateFixedTerms(system::System &to_update,
                        system::EnergyGroup group,
                        quadrature::QuadraturePointIndex index) override;
  void UpdateFissionSource(system::System &to_update,
                           system::EnergyGroup group,
                           quadrature::QuadraturePointIndex index) override;
  void UpdateScatteringSource(system::System &to_update,
                              system::EnergyGroup group,
                              quadrature::QuadraturePointIndex index) override;

  EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map() const {
    return angular_solution_ptr_map_; }
  std::unordered_set<Boundary> reflective_boundaries() const {
    return reflective_boundaries_; }
  SAAFFormulationType* formulation_ptr() const {return formulation_ptr_.get();};
  StamperType* stamper_ptr() const {return stamper_ptr_.get();};
  QuadratureSetType* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get();};
 private:
  bool IsOnReflectiveBoundary(const domain::CellPtr<dim>& cell_ptr,
                              const domain::FaceIndex face_index) const {
    return reflective_boundaries_.count(
        static_cast<const Boundary>(
            cell_ptr->face(face_index.get())->boundary_id()));
  }
  std::unique_ptr<SAAFFormulationType> formulation_ptr_;
  std::unique_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map_;
  std::unordered_set<Boundary> reflective_boundaries_ = {};
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_SAAF_UPDATER_H_

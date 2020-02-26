#ifndef BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_
#define BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_

#include "domain/definition_i.h"
#include "formulation/angular_stamper_i.h"

#include "formulation/angular/self_adjoint_angular_flux_i.h"

#include <memory>

namespace bart {

namespace formulation {

template <int dim>
class CFEM_SAAF_Stamper : public AngularStamperI<dim> {
 public:
  static constexpr int dimension = dim;
  using DomainDefinitionType = domain::DefinitionI<dim>;
  using SAAFFormulationType = typename
      formulation::angular::SelfAdjointAngularFluxI<dim>;

  CFEM_SAAF_Stamper(std::unique_ptr<SAAFFormulationType> saaf_ptr,
                    std::shared_ptr<DomainDefinitionType> defintion_ptr);

  void StampBoundaryBilinearTerm(
      system::MPISparseMatrix &to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;

  void StampCollisionTerm(system::MPISparseMatrix& to_stamp,
                          const system::EnergyGroup group_number) override;

  void StampFissionSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const double k_eff,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments) override;

  void StampFixedSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;

  void StampScatteringSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments) override;

  void StampStreamingTerm(
      system::MPISparseMatrix& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) override;

  DomainDefinitionType* definition_ptr() const {return definition_ptr_.get();}
  SAAFFormulationType* formulation_ptr() const {return formulation_ptr_.get();}

 private:

  // Dependencies
  std::unique_ptr<SAAFFormulationType> formulation_ptr_ = nullptr;
  std::shared_ptr<DomainDefinitionType> definition_ptr_ = nullptr;

  std::vector<domain::CellPtr<dim>> cells_;

  // Private functions
  void StampMatrix(
      system::MPISparseMatrix& to_stamp,
      std::function<void(formulation::FullMatrix&,
                         const domain::CellPtr<dim>&)> stamping_function);
  void StampVector(
      system::MPIVector& to_stamp,
      std::function<void(formulation::Vector&,
                         const domain::CellPtr<dim>&)> stamping_function);
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_CFEM_SAAF_STAMPER_H_

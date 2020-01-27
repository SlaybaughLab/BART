#ifndef BART_SRC_FORMULATION_ANGULAR_STAMPER_I_H_
#define BART_SRC_FORMULATION_ANGULAR_STAMPER_I_H_

#include <memory>

#include "formulation/stamper_i.h"
#include "quadrature/quadrature_point_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/system_types.h"

namespace bart {

namespace formulation {

template <int dim>
class AngularStamperI : public StamperI {
 public:
  static constexpr int dimension = dim;
  virtual ~AngularStamperI() = default;
  virtual void StampCollisionTerm(system::MPISparseMatrix& to_stamp,
                                  const system::EnergyGroup group_number) = 0;
  virtual void StampFissionSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
  const system::EnergyGroup group_number,
  const double k_eff,
  const system::moments::MomentVector &in_group_moment,
  const system::moments::MomentsMap &group_moments) = 0;

  virtual void StampFixedSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) = 0;

  virtual void StampScatteringSourceTerm(
      system::MPIVector& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number,
      const system::moments::MomentVector &in_group_moment,
      const system::moments::MomentsMap &group_moments) = 0;

  virtual void StampStreamingTerm(
      system::MPISparseMatrix& to_stamp,
      const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
      const system::EnergyGroup group_number) = 0;
};

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_STAMPER_I_H_

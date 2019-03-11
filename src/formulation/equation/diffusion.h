#ifndef BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
#define BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_

#include "formulation/types.h"
#include "formulation/equation/transport_scalar.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Diffusion : public TransportScalar<dim> {
 public:
  using MaterialID = int;
  using typename TransportScalar<dim>::CellPtr;
  using typename TransportScalar<dim>::Matrix;
  using typename TransportScalar<dim>::GroupNumber;

  explicit Diffusion(const DiscretizationType discretization)
      : TransportScalar<dim>(discretization) {}

  virtual ~Diffusion() = default;

  void Precalculate() override;

  void FillCellBilinearTerm(Matrix& to_fill,
                            const CellPtr &cell_ptr,
                            const GroupNumber group) const override;



 protected:
  using TransportScalar<dim>::finite_element_;
  using TransportScalar<dim>::cross_sections_;
  using TransportScalar<dim>::cell_degrees_of_freedom_;
  using TransportScalar<dim>::cell_quadrature_points_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
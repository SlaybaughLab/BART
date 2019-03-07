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
  using typename TransportScalar<dim>::CellPtr;
  using typename TransportScalar<dim>::Matrix;
  using typename TransportScalar<dim>::GroupNumber;

  explicit Diffusion(const DiscretizationType discretization)
      : TransportScalar<dim>(discretization) {}

  virtual ~Diffusion() = default;

  void FillCellBilinearTerm(Matrix& to_fill,
                            const CellPtr &cell_ptr,
                            const GroupNumber group) override;
 protected:
  using TransportScalar<dim>::finite_element_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
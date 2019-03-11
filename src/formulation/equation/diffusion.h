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
  using typename TransportScalar<dim>::FaceNumber;

  explicit Diffusion(const DiscretizationType discretization)
      : TransportScalar<dim>(discretization) {}

  virtual ~Diffusion() = default;

  void Precalculate(const CellPtr &cell_ptr) override;

  void FillCellBilinearTerm(Matrix& to_fill,
                            const CellPtr &cell_ptr,
                            const GroupNumber group) const override;

  void FillBoundaryBilinearTerm(Matrix &to_fill,
                                const CellPtr &cell_ptr,
                                const FaceNumber face_number) const override;

  void FillCellLinearScatteringTerm(Matrix &to_fill,
                                    const CellPtr &cell_ptr,
                                    const GroupNumber group) const override;


 protected:

  std::vector<Matrix> shape_squared_;
  std::vector<Matrix> gradient_squared_;

  using TransportScalar<dim>::scalar_fluxes_;
  using TransportScalar<dim>::SetCell;
  using TransportScalar<dim>::SetFace;
  using TransportScalar<dim>::finite_element_;
  using TransportScalar<dim>::cross_sections_;
  using TransportScalar<dim>::cell_degrees_of_freedom_;
  using TransportScalar<dim>::cell_quadrature_points_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
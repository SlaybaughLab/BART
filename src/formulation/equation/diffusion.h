#ifndef BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
#define BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_

#include "problem/parameter_types.h"
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
  using typename TransportScalar<dim>::Vector;
  using typename TransportScalar<dim>::GroupNumber;
  using typename TransportScalar<dim>::FaceNumber;

  Diffusion(const DiscretizationType discretization,
            problem::ProblemType problem_type,
            std::shared_ptr<double> k_effective)
      : TransportScalar<dim>(discretization),
        problem_type_(problem_type),
        k_effective_(k_effective) {}

  virtual ~Diffusion() = default;

  void Precalculate(const CellPtr &cell_ptr) override;

  void FillCellFixedBilinear(Matrix &to_fill,
                             const CellPtr &cell_ptr,
                             const GroupNumber group) const override;

  void FillBoundaryFixedBilinear(Matrix &to_fill,
                                 const CellPtr &cell_ptr,
                                 const GroupNumber group,
                                 const FaceNumber face_number,
                                 const BoundaryType boundary_type) const override;

  void FillCellFixedLinear(Vector &rhs_to_fill,
                           const CellPtr &cell_ptr,
                           const GroupNumber group) const override;

  void FillCellVariableLinear(Vector& rhs_to_fill,
                              const CellPtr &cell_ptr,
                              const GroupNumber group) const override;

 protected:
  problem::ProblemType problem_type_;
  std::shared_ptr<double> k_effective_;
  std::vector<Matrix> shape_squared_;
  std::vector<Matrix> gradient_squared_;

  using TransportScalar<dim>::scalar_fluxes_;
  using TransportScalar<dim>::SetCell;
  using TransportScalar<dim>::SetFace;
  using TransportScalar<dim>::finite_element_;
  using TransportScalar<dim>::cross_sections_;
  using TransportScalar<dim>::cell_degrees_of_freedom_;
  using TransportScalar<dim>::cell_quadrature_points_;
  using TransportScalar<dim>::face_quadrature_points_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
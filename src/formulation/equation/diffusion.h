#ifndef BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
#define BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_

#include "problem/parameter_types.h"
#include "formulation/types.h"
#include "formulation/equation/scalar_fixed_bilinear.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Diffusion : public ScalarFixedBilinear<dim> {
 public:
  using MaterialID = int;
  using typename ScalarFixedBilinear<dim>::CellPtr;
  using typename ScalarFixedBilinear<dim>::Matrix;
  using typename ScalarFixedBilinear<dim>::Vector;
  using typename ScalarFixedBilinear<dim>::GroupNumber;
  using typename ScalarFixedBilinear<dim>::FaceNumber;

  Diffusion(const DiscretizationType discretization,
            problem::ProblemType problem_type,
            std::shared_ptr<double> k_effective)
      : ScalarFixedBilinear<dim>(discretization),
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

  using ScalarFixedBilinear<dim>::scalar_fluxes_;
  using ScalarFixedBilinear<dim>::SetCell;
  using ScalarFixedBilinear<dim>::SetFace;
  using ScalarFixedBilinear<dim>::finite_element_;
  using ScalarFixedBilinear<dim>::cross_sections_;
  using ScalarFixedBilinear<dim>::cell_degrees_of_freedom_;
  using ScalarFixedBilinear<dim>::cell_quadrature_points_;
  using ScalarFixedBilinear<dim>::face_quadrature_points_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_DIFFUSION_H_
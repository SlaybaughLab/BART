#ifndef BART_SRC_FORMULATION_EQUATION_SCALAR_FIXED_BILINEAR_H_
#define BART_SRC_FORMULATION_EQUATION_SCALAR_FIXED_BILINEAR_H_

#include "data/forward_declarations.h"
#include "formulation/types.h"
#include "formulation/equation/transport.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class ScalarFixedBilinear : public Transport<dim> {
 public:
  using typename Transport<dim>::CellPtr;
  using Matrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using GroupNumber = int;
  using FaceNumber = int;

  virtual ~ScalarFixedBilinear() = default;

  ScalarFixedBilinear(const DiscretizationType discretization)
      : Transport<dim>(EquationType::kScalar, discretization) {};

  /*! \brief Pre-calculates any part of the formulation that is not cell and group
   * dependent.
   *
   * \param[in] cell_ptr any cell in the current mesh. This is required to
   * initialize the finite element object for the gradient and shape values
   * (cell differences are determined by Jacobians and materials). This does not
   * need to be called for every cell, just once at the beginning of assembly.
   */
  virtual void Precalculate(const CellPtr &cell_ptr) = 0;

  /*! \brief Fills a cell matrix with the integrated fixed bilinear term for a
   * cell.
   *
   * These are the bilinear terms that will not change between iterations.
   *
   * \param[in,out] to_fill local cell matrix to fill with the integrated term.
   * \param[in] cell_ptr the cell to use to generate the local cell matrix.
   * \param[in] group energy group number
   */
  virtual void FillCellFixedBilinear(Matrix& to_fill,
                                     const CellPtr &cell_ptr,
                                     const GroupNumber group) const = 0;

  /*! \brief Fills a cell matrix with the integrated fixed boundary bilinear
   * term for a boundary cell.
   *
   * \param[in, out] to_fill local cell matrix to fill with the integrated term.
   * \param[in] cell_ptr the boundary cell used to generate the local cell matrix.
   * \param[in] group energy group number
   * \param[in] face_number the face number for the boundary.
   * \param[in] boundary_type the type of boundary.
   */
  virtual void FillBoundaryFixedBilinear(
      Matrix &to_fill,
      const CellPtr &cell_ptr,
      const GroupNumber group,
      const FaceNumber face_number,
      const BoundaryType boundary_type) const = 0;

  /*! \brief Fills a cell vector with the integrated fixed linear terms.
   *
   * \param[in, out] rhs_to_fill vector to fill with the integrated term.
   * \param[in] cell_ptr the cell used to generate the local cell matrix.
   * \param[in] group energy group number
   */
  virtual void FillCellFixedLinear(Vector& rhs_to_fill,
                                   const CellPtr &cell_ptr,
                                   const GroupNumber group) const = 0;

  /*! \brief Fills a cell vector with the integrated variable linear term.
   *
   * \param[in, out] rhs_to_fill vector to fill with the integrated term.
   * \param[in] cell_ptr the cell used to generate the local cell matrix.
   * \param[in] group energy group number
   */
  virtual void FillCellVariableLinear(
      Vector& rhs_to_fill,
      const CellPtr &cell_ptr,
      const GroupNumber group,
      const data::ScalarFluxPtrs &scalar_flux) const = 0;

  ScalarEquations scalar_equation() const { return scalar_equation_; };
 protected:
  ScalarEquations scalar_equation_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_
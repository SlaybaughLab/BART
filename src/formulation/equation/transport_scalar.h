#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_

#include "formulation/types.h"
#include "formulation/equation/transport.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class TransportScalar : public Transport<dim> {
 public:
  using typename Transport<dim>::CellPtr;
  using Matrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using GroupNumber = int;
  using FaceNumber = int;

  virtual ~TransportScalar() = default;

  TransportScalar(const DiscretizationType discretization)
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

  /*! \brief Fills a cell matrix with the integrated bilinear term for a cell.
   *
   * \param[in,out] to_fill local cell matrix to fill with the integrated term.
   * \param[in] cell_ptr the cell to use to generate the local cell matrix.
   * \param[in] group energy group number
   */
  virtual void FillCellBilinearTerm(Matrix& to_fill,
                                    const CellPtr &cell_ptr,
                                    const GroupNumber group) const = 0;
  /*! \brief Fills a cell matrix with the integrated boundary term for a
   * boundary cell.
   *
   * \param[in, out] to_fill local cell matrix to fill with the integrated term.
   * \param[in] cell_ptr the boundary cell used to generate the local cell matrix.
   * \param[in] group energy group number
   * \param[in] face_number the face number for the boundary.
   * \param boundary_type the type of boundary.
   */
  virtual void FillBoundaryBilinearTerm(
      Matrix &to_fill,
      const CellPtr &cell_ptr,
      const GroupNumber group,
      const FaceNumber face_number,
      const BoundaryType boundary_type) const = 0;

  virtual void FillCellLinearScatteringTerm(Vector& rhs_to_fill,
                                            const CellPtr &cell_ptr,
                                            const GroupNumber group) const = 0;

  ScalarEquations scalar_equation() const { return scalar_equation_; };
 protected:
  ScalarEquations scalar_equation_;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_SCALAR_H_
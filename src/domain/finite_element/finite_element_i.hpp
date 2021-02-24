#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_I_H
#define BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

#include <deal.II/fe/fe_values.h>

#include "domain/domain_types.h"
#include "utility/has_description.h"

//! Classes that provide a finite-element basis
namespace bart::domain::finite_element {

/*! \brief Interface for a finite element object based on the dealii library.
 *
 * This object is responsible for returning finite element data such as: basis
 * function values at quadrature points, basis function values on faces, numbers
 * of quadrature points and dofs per cell.
 *
 */
template <int dim>
class FiniteElementI : public utility::HasDescription {
 public:
  using CellPtr = domain::CellPtr<dim>;
  using DealiiVector = dealii::Vector<double>;
  using Tensor = dealii::Tensor<1, dim>;
  virtual ~FiniteElementI() = default;

  // Basic FE properties
  /*! \brief Gets polynomial degree */
  virtual auto polynomial_degree() const -> int = 0;
  /*! \brief Gets degrees of freedom per cell */
  virtual auto dofs_per_cell() const -> int = 0;
  /*! \brief Gets number of quadrature points per cell */
  virtual auto n_cell_quad_pts() const -> int = 0;
  /*! \brief Gets number of quadrature points per face */
  virtual auto n_face_quad_pts() const -> int = 0;

  // Various methods to access the underlying finite element object
  /*! \brief Sets the cell for finite element values.
   *
   * Does nothing if the cell is already set to the passed value.
   *
   * \param to_set cell to set
   * \return bool indicating if the cell was changed.
   */
  virtual auto SetCell(const CellPtr &to_set) -> bool = 0;

  /*! \brief Sets the face and cell.
   *
   * \param to_set cell to set
   * \param face face number to set
   * \return indicating if the cell was changed.
   */
  virtual auto SetFace(const CellPtr &to_set, const domain::FaceIndex face) -> bool = 0;

  /*! \brief Get the value of shape functions.
   *
   * Returns the value of the \f$i\f$-th shape function, \f$\varphi_i\f$, that
   * corresponds to the \f$i\f$-th cell degree of freedom, evaluated at cell
   * quadrature point, \f$q\f$.
   *
   * \param cell_degree_of_freedom value of \f$i\f$.
   * \param cell_quadrature_point cell quadrature point \f$q\f$ to evaluate the
   * shape function.
   * \return double corresponding to the value.
   */
  virtual auto ShapeValue(const int cell_degree_of_freedom, const int cell_quadrature_point) const -> double = 0;

  /*! \brief Get the value of face shape functions.
   *
   * Returns the value of the \f$i\f$-th shape function, \f$\varphi_i\f$, that
   * corresponds to the \f$i\f$-th cell degree of freedom, evaluated at face
   * quadrature point, \f$q\f$.
   *
   * \param cell_degree_of_freedom value of \f$i\f$.
   * \param face_quadrature_point face quadrature point \f$q\f$ to evaluate the
   * shape function.
   * \return double corresponding to the value.
   */
  virtual auto FaceShapeValue(const int cell_degree_of_freedom, const int face_quadrature_point) const -> double = 0;

  /*! \brief Get the value of gradient of the shape function.
   *
   * Returns the value of the gradient of the \f$i\f$-th shape function,
   * \f$\nabla\varphi_i\f$, that corresponds to the \f$i\f$-th cell degree of
   * freedom, evaluated at cell quadrature point, \f$q\f$.
   *
   * \param cell_degree_of_freedom value of \f$i\f$.
   * \param cell_quadrature_point cell quadrature point \f$q\f$ to evaluate the
   * shape function.
   * \return dealii::Tensor with the value of the gradient.
   */
  virtual auto ShapeGradient(const int cell_degree_of_freedom, const int cell_quadrature_point) const -> Tensor = 0;

  /*! \brief Get the value of the Jacobian for the current cell.
   *
   * Returns the Jacobian evaluated at the specified quadrature point, needed
   * for all integrations.
   *
   * \param cell_quadrature_point cell quadrature point \f$q\f$ to evaluate the
   * Jacobian function.
   * \return double corresponding to the value.
   */
  virtual auto Jacobian(const int cell_quadrature_point) const -> double = 0;

  /*! \brief Get the value of the Jacobian for the current face.
 *
 * Returns the Jacobian evaluated at the specified face quadrature point, needed
 * for all integrations.
 *
 * \param face_quadrature_point face quadrature point \f$q\f$ to evaluate the
 * Jacobian function.
 * \return double corresponding to the value.
 */
  virtual auto FaceJacobian(const int face_quadrature_point) const -> double = 0;

  /*! \brief Get the value of the face normal for the current face.
   *
   * Returns the normal vector for the specified face.
   *
   */
  virtual auto FaceNormal() const -> Tensor = 0;

  /*! \brief Get the value of a vector at the interior cell quadrature points.
   *
   * \param vector_at_dofs vector value at dofs
   * \return a vector holding the value of the vector at each quadrature point.
   */
  virtual auto ValueAtQuadrature(const DealiiVector& vector_at_dofs) const -> std::vector<double> = 0;

  /*! \brief Get the value of an MPI Vector at the cell face quadrature points.
   *
   * @param mpi_vector mpi vector to get the face values of.
   * @return a vector holding the value of the mpi vector at each face quadrature point.
   */
   virtual auto ValueAtFaceQuadrature(const DealiiVector& values_at_dofs) const -> std::vector<double> = 0;

  // DealII Finite element object access. These methods access the underlying
  // finite element objects.
  /*! \brief Gets pointer to the underlying finite element object */
  virtual auto finite_element() -> dealii::FiniteElement<dim, dim>* = 0;
  /*! \brief Gets pointer to the underlying finite element values object */
  virtual auto values() -> dealii::FEValues<dim>* = 0;
  /*! \brief Gets pointer to underlying finite element faces values object */
  virtual auto face_values() -> dealii::FEFaceValues<dim>* = 0;
  /*! \brief Gets pointer to second face values object.
   * This is often used if you need to get face values for a neighbor cell. */
  virtual auto neighbor_face_values() -> dealii::FEFaceValues<dim>* = 0;
};

} // namespace bart::domain::finite_element

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

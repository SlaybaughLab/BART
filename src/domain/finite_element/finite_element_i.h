#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_I_H
#define BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

#include <deal.II/fe/fe_values.h>

#include "domain/domain_types.h"
#include "system/moments/spherical_harmonic_types.h"

/*! \brief Interface for a finite element object based on the dealii library.
 *
 * This object is responsible for returning finite element data such as: basis
 * function values at quadrature points, basis function values on faces, numbers
 * of quadrature points and dofs per cell.
 * 
 */

namespace bart {

namespace domain {

namespace finite_element {

template <int dim>
class FiniteElementI {
 public:
  virtual ~FiniteElementI() = default;

  // Basic FE properties
  /*! \brief Gets polynomial degree */
  virtual int polynomial_degree() const = 0;
  /*! \brief Gets degrees of freedom per cell */
  virtual int dofs_per_cell() const = 0;
  /*! \brief Gets number of quadrature points per cell */
  virtual int n_cell_quad_pts() const = 0;
  /*! \brief Gets number of quadrature points per face */
  virtual int n_face_quad_pts() const = 0;

  // Various methods to access the underlying finite element object
  /*! \brief Sets the cell for finite element values.
   *
   * Does nothing if the cell is already set to the passed value.
   *
   * \param to_set cell to set
   * \return bool indicating if the cell was changed.
   */
  virtual bool SetCell(const domain::CellPtr<dim> &to_set) = 0;

  /*! \brief Sets the face and cell.
   *
   * \param to_set cell to set
   * \param face face number to set
   * \return indicating if the cell was changed.
   */
  virtual bool SetFace(const domain::CellPtr<dim> &to_set,
                       const domain::FaceIndex face) = 0;

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
  virtual double ShapeValue(const int cell_degree_of_freedom,
                            const int cell_quadrature_point) const = 0;

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
  virtual double FaceShapeValue(const int cell_degree_of_freedom,
                                const int face_quadrature_point) const = 0;

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
  virtual dealii::Tensor<1, dim> ShapeGradient(
      const int cell_degree_of_freedom,
      const int cell_quadrature_point) const = 0;

  /*! \brief Get the value of the Jacobian for the current cell.
   *
   * Returns the Jacobian evaluated at the specified quadrature point, needed
   * for all integrations.
   *
   * \param cell_quadrature_point cell quadrature point \f$q\f$ to evaluate the
   * Jacobian function.
   * \return double corresponding to the value.
   */
  virtual double Jacobian(const int cell_quadrature_point) const = 0;

  /*! \brief Get the value of the Jacobian for the current face.
 *
 * Returns the Jacobian evaluated at the specified face quadrature point, needed
 * for all integrations.
 *
 * \param face_quadrature_point face quadrature point \f$q\f$ to evaluate the
 * Jacobian function.
 * \return double corresponding to the value.
 */
  virtual double FaceJacobian(const int face_quadrature_point) const = 0;

  /*! \brief Get the value of the face normal for the current face.
   *
   * Returns the normal vector for the specified face.
   *
   */
   virtual dealii::Tensor<1, dim> FaceNormal() const = 0;

  /*! \brief Get the value of a flux moment at the interior cell quadrature points.
   *
   * \param moment flux moment to get the value of.
   * \return a vector holding the value of the moment at each quadrature point.
   */
  virtual std::vector<double> ValueAtQuadrature(
      const system::moments::MomentVector moment) const = 0;

  // DealII Finite element object access. These methods access the underlying
  // finite element objects.
  /*! \brief Gets pointer to the underlying finite element object */
  virtual dealii::FiniteElement<dim, dim> *finite_element() = 0;
  /*! \brief Gets pointer to the underlying finite element values object */
  virtual dealii::FEValues<dim> *values() = 0;
  /*! \brief Gets pointer to underlying finite element faces values object */
  virtual dealii::FEFaceValues<dim> *face_values() = 0;
  /*! \brief Gets pointer to second face values object.
   * This is often used if you need to get face values for a neighbor cell. */
  virtual dealii::FEFaceValues<dim> *neighbor_face_values() = 0;
};

} // namespace finite_element

} // namespace domain

} // namespace bart 

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

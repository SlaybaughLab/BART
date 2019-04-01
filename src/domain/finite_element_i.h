#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_I_H
#define BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

#include <deal.II/fe/fe_values.h>

/*! \brief Interface for a finite element object based on the dealii library.
 *
 * This object is responsible for returning finite element data such as: basis
 * function values at quadrature points, basis function values on faces, numbers
 * of quadrature points and dofs per cell.
 * 
 */

namespace bart {

namespace domain {

template <int dim>
class FiniteElementI {
 public:
  virtual ~FiniteElementI() = default;
  using FaceNumber = int;
  using CellPtr = typename dealii::DoFHandler<dim>::active_cell_iterator;

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
  virtual bool SetCell(const CellPtr &to_set) = 0;

  /*! \brief Sets the face and cell.
   *
   * \param to_set cell to set
   * \param face face number to set
   * \return indicating if the cell was changed.
   */
  virtual bool SetFace(const CellPtr &to_set,
                       const FaceNumber face) = 0;

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

} // namespace domain

} // namespace bart 

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_I_H

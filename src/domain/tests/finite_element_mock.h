#ifndef BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_H
#define BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_H

#include "../finite_element_i.h"

#include <deal.II/fe/fe_values.h>

#include "../../test_helpers/gmock_wrapper.h"

namespace bart {

namespace domain {

template <int dim>
class FiniteElementMock : public FiniteElementI<dim> {
 public:

  using typename FiniteElementI<dim>::CellPtr;
  using typename FiniteElementI<dim>::FaceNumber;

  MOCK_CONST_METHOD0_T(polynomial_degree, int());

  MOCK_CONST_METHOD0_T(dofs_per_cell, int());

  MOCK_CONST_METHOD0_T(n_cell_quad_pts, int());

  MOCK_CONST_METHOD0_T(n_face_quad_pts, int());

  MOCK_METHOD1_T(SetCell, bool(const CellPtr &));

  MOCK_METHOD2_T(SetFace, bool(const CellPtr &to_set, const FaceNumber face));

  MOCK_CONST_METHOD2_T(ShapeValue, double(const int, const int));

  MOCK_CONST_METHOD2_T(ShapeGradient, dealii::Tensor<1, dim>(const int, const int));

  MOCK_CONST_METHOD1_T(Jacobian, double(const int));

  MOCK_METHOD0_T(finite_element, dealii::FiniteElement<dim, dim>*());

  MOCK_METHOD0_T(values, dealii::FEValues<dim>*());

  MOCK_METHOD0_T(face_values, dealii::FEFaceValues<dim>*());

  MOCK_METHOD0_T(neighbor_face_values, dealii::FEFaceValues<dim>*());
};

} // namespace domain

} // namespace bart 

#endif //BART_SRC_DOMAIN_FINITE_ELEMENT_MOCK_H

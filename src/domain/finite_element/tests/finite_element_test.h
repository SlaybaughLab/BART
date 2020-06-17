#ifndef BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
#define BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_

#include "domain/finite_element/finite_element.h"

#include "test_helpers/test_assertions.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace bart {

namespace domain {

namespace finite_element {

namespace testing {

template <int dim>
class FiniteElementBaseClassTest : public ::testing::Test {
 public:
  FiniteElementBaseClassTest()
      : dof_handler_(triangulation_) {}
  virtual ~FiniteElementBaseClassTest() = default;
 protected:
  dealii::Triangulation<dim> triangulation_;
  dealii::DoFHandler<dim> dof_handler_;

  void TestSetCell(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestSetCellAndFace(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestValueAtQuadrature(domain::finite_element::FiniteElement<dim>* test_fe);
  void TestValueAtFaceQuadrature(domain::finite_element::FiniteElement<dim>* test_fe);
  void SetUp() override {
    dealii::GridGenerator::hyper_cube(triangulation_, -1, 1);
    triangulation_.refine_global(2);
  }
};

template <int dim>
void FiniteElementBaseClassTest<dim>::TestSetCell(
    FiniteElement<dim> *test_fe) {

  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();
  auto cell_id = cell->id();

  EXPECT_NO_THROW(test_fe->SetCell(cell));

  test_fe->values()->reinit(cell);

  EXPECT_FALSE(test_fe->SetCell(cell)); // Shouldn't change anything
  EXPECT_EQ(cell_id, test_fe->values()->get_cell()->id()); // Cell didn't change

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();

  EXPECT_TRUE(test_fe->SetCell(next_cell));
  // Check changed
  EXPECT_NE(cell_id, test_fe->values()->get_cell()->id());
  EXPECT_EQ(next_cell_id, test_fe->values()->get_cell()->id());

}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestSetCellAndFace(
    FiniteElement<dim> *test_fe) {
  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();
  auto cell_id = cell->id();
  int face = 0;
  int face_index = cell->face_index(face);

  test_fe->SetFace(cell, domain::FaceIndex(face));

  test_fe->face_values()->reinit(cell, face);

  EXPECT_FALSE(test_fe->SetFace(cell, domain::FaceIndex(face)));
  EXPECT_EQ(cell_id, test_fe->face_values()->get_cell()->id());
  EXPECT_EQ(face_index, test_fe->face_values()->get_face_index());

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();
  int next_face = face + 1;
  int next_face_index = next_cell->face_index(next_face);

  EXPECT_TRUE(test_fe->SetFace(next_cell, domain::FaceIndex(next_face)));
  EXPECT_EQ(next_cell_id, test_fe->face_values()->get_cell()->id());
  EXPECT_NE(face_index, test_fe->face_values()->get_face_index());
  EXPECT_EQ(next_face_index, test_fe->face_values()->get_face_index());
}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestValueAtQuadrature(
    FiniteElement<dim> *test_fe) {

  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();

  EXPECT_NO_THROW(test_fe->SetCell(cell));

  int n_dofs = dof_handler_.n_dofs();

  std::vector<double> moment_values(n_dofs, 0.5);
  system::moments::MomentVector test_moment(moment_values.begin(), moment_values.end());

  std::vector<double> expected_vector(test_fe->n_cell_quad_pts(), 0.5);

  auto result_vector = test_fe->ValueAtQuadrature(test_moment);

  EXPECT_TRUE(bart::test_helpers::CompareVector(expected_vector, result_vector));
}

template <int dim>
void FiniteElementBaseClassTest<dim>::TestValueAtFaceQuadrature(
    FiniteElement<dim> *test_fe) {

  dof_handler_.distribute_dofs(*test_fe->finite_element());

  auto cell = dof_handler_.begin_active();

  EXPECT_NO_THROW(test_fe->SetFace(cell, domain::FaceIndex(0)));

  int n_dofs = dof_handler_.n_dofs();

  system::MPIVector test_vector(MPI_COMM_WORLD, n_dofs, n_dofs);
  test_vector = 0.5;
  test_vector.compress(dealii::VectorOperation::insert);

  std::vector<double> expected_vector(test_fe->n_face_quad_pts(), 0.5);

  auto result_vector = test_fe->ValueAtFaceQuadrature(test_vector);

  EXPECT_TRUE(bart::test_helpers::CompareVector(expected_vector, result_vector));
}

template <int dim>
class DomainFiniteElementBaseDomainTest :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<dim> {
 public:
  void SetUp() override;
  void TestValueAtFaceQuadrature(FiniteElement<dim> *test_fe);
};

template <int dim>
void DomainFiniteElementBaseDomainTest<dim>::SetUp() {
  this->SetUpDealii();
}

template <int dim>
void DomainFiniteElementBaseDomainTest<dim>::TestValueAtFaceQuadrature(
    FiniteElement<dim> *test_fe) {
  EXPECT_TRUE(false);
}



} // namespace testing

} // namespace finite_element

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
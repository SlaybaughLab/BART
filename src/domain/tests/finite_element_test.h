#ifndef BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
#define BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_

#include "domain/finite_element.h"

#include "test_helpers/test_assertions.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace domain {

namespace testing {

template <int dim>
class FiniteElementBaseClassTest : public ::testing::Test {
 public:
  FiniteElementBaseClassTest()
      : dof_handler_(triangulation_) {}
 protected:
  dealii::Triangulation<dim> triangulation_;
  dealii::DoFHandler<dim> dof_handler_;

  void TestSetCell(domain::FiniteElement<dim>* test_fe);
  void TestSetCellAndFace(domain::FiniteElement<dim>* test_fe);
  void TestValueAtQuadrature(domain::FiniteElement<dim>* test_fe);
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

  test_fe->SetFace(cell, face);

  test_fe->face_values()->reinit(cell, face);

  EXPECT_FALSE(test_fe->SetFace(cell, face));
  EXPECT_EQ(cell_id, test_fe->face_values()->get_cell()->id());
  EXPECT_EQ(face_index, test_fe->face_values()->get_face_index());

  auto next_cell = cell;
  ++next_cell;
  auto next_cell_id = next_cell->id();
  int next_face = face + 1;
  int next_face_index = next_cell->face_index(next_face);

  EXPECT_TRUE(test_fe->SetFace(next_cell, next_face));
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

  std::vector<double> expected_vector(test_fe->dofs_per_cell(), 0.5);

  auto result_vector = test_fe->ValueAtQuadrature(test_moment);

  EXPECT_TRUE(bart::testing::CompareVector(expected_vector, result_vector));
}

} // namespace testing

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
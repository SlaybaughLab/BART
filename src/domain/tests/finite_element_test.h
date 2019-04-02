#ifndef BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
#define BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_

#include "domain/finite_element.h"

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

};

template <int dim>
void FiniteElementBaseClassTest<dim>::TestSetCell(
    FiniteElement<dim> *test_fe) {
  dof_handler_.distribute_dofs(test_fe->finite_element());

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

} // namespace testing

} // namespace domain

} // namespace bart

#endif // BART_SRC_DOMAIN_TESTS_FINITE_ELEMENT_TEST_H_
#include "data/system_scalar_fluxes.h"

#include <gtest/gtest.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include "test_helpers/gmock_wrapper.h"

class DataSystemScalarFluxesTest : public ::testing::Test {};

TEST_F(DataSystemScalarFluxesTest, BuilderTest){
  dealii::Triangulation<2> triangulation;
  dealii::FE_Q<2> fe(1);
  dealii::GridGenerator::hyper_cube (triangulation, -1, 1);
  dealii::DoFHandler<2> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  int n_groups = 3;
  int n_dofs = dof_handler.n_dofs();
  auto locally_owned_dofs = dof_handler.locally_owned_dofs();

  auto flux_ptrs_ptr = bart::data::BuildSystemScalarFluxes(
      n_groups, locally_owned_dofs);

  EXPECT_EQ(flux_ptrs_ptr->size(), n_groups);
  for (const auto &flux_pair : *flux_ptrs_ptr) {
    EXPECT_EQ((flux_pair.second)->size(), n_dofs);
  }
}
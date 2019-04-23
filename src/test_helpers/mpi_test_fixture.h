#ifndef BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_
#define BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_

#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace testing {

template <int dim>
 class MPI_TestFixture : public ::testing::Test {
  protected:
   MPI_TestFixture()
   : triangulation_(MPI_COMM_WORLD,
                    typename dealii::Triangulation<2>::MeshSmoothing(
                        dealii::Triangulation<2>::smoothing_on_refinement |
                            dealii::Triangulation<2>::smoothing_on_coarsening)),
     dof_handler_(triangulation_),
     fe_(1) {};

   void SetUpDealii();

   using Cell = typename dealii::DoFHandler<dim>::active_cell_iterator;
   dealii::ConstraintMatrix constraint_matrix_;
   dealii::parallel::distributed::Triangulation<dim> triangulation_;
   dealii::DoFHandler<dim> dof_handler_;
   dealii::FE_Q<dim> fe_;
   dealii::IndexSet locally_relevant_dofs;
   dealii::IndexSet locally_owned_dofs_;
   std::vector<Cell> cells_;

   dealii::PETScWrappers::MPI::SparseMatrix matrix_1, matrix_2, matrix_3;
};

template <int dim>
inline void MPI_TestFixture<dim>::SetUpDealii() {
  dealii::GridGenerator::hyper_cube(triangulation_, 0, 1);
  dof_handler_.distribute_dofs(fe_);
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_,
                                                  locally_relevant_dofs);
  locally_owned_dofs_ = dof_handler_.locally_owned_dofs();

  for (auto cell = dof_handler_.begin_active(); cell != dof_handler_.end(); ++ cell) {
    if (cell->is_locally_owned())
      cells_.push_back(cell);
  }

  constraint_matrix_.clear();
  constraint_matrix_.reinit(locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();

  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp,
                                          constraint_matrix_, false);
  dealii::SparsityTools::distribute_sparsity_pattern(
      dsp,
      dof_handler_.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs);
  constraint_matrix_.condense(dsp);

  matrix_1.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp, MPI_COMM_WORLD);
  matrix_2.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp, MPI_COMM_WORLD);
  matrix_3.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp, MPI_COMM_WORLD);
  }

} // namespace testing

} // namespace bart

#endif // BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_

#include <deal.II/base/mpi.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "formulation/assembler.h"
#include "formulation/scalar/cfem_diffusion.h"

#include "test_helpers/gmock_wrapper.h"

class FormulationAssemblerCFEMDiffusionTest : public ::testing::Test {

};

TEST_F(FormulationAssemblerCFEMDiffusionTest, MPIDummy) {
  int this_mpi_process(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
  int n_mpi_processes(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD));
  dealii::ConditionalOStream pcout(std::cout, this_mpi_process == 0);
  dealii::ConstraintMatrix constraint_matrix;
  typedef std::vector<typename dealii::DoFHandler<2>::active_cell_iterator> CellRange;
  CellRange cells;

  dealii::parallel::distributed::Triangulation<2> triangulation (
      MPI_COMM_WORLD,
      typename dealii::Triangulation<2>::MeshSmoothing(
          dealii::Triangulation<2>::smoothing_on_refinement |
          dealii::Triangulation<2>::smoothing_on_coarsening));
  dealii::GridGenerator::hyper_cube(triangulation, 0, 1);
  triangulation.refine_global(2);

  dealii::DoFHandler<2> dof_handler(triangulation);
  dealii::FE_Q<2> fe(1);
  dof_handler.distribute_dofs(fe);

  dealii::IndexSet locally_relevant_dofs;
  auto locally_owned_dofs = dof_handler.locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler,
                                                  locally_relevant_dofs);

  constraint_matrix.clear();
  constraint_matrix.reinit(locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler,
                                                  constraint_matrix);
  constraint_matrix.close();

  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);

  dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp,
                                          constraint_matrix, false);
  dealii::SparsityTools::distribute_sparsity_pattern(dsp,
      dof_handler.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs);
  constraint_matrix.condense(dsp);

  dealii::PETScWrappers::MPI::SparseMatrix system_matrix;
  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  for (auto cell = dof_handler.begin_active(); cell != dof_handler.end(); ++ cell) {
    if (cell->is_locally_owned())
      cells.push_back(cell);
  }

  dealii::FullMatrix<double> cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
  for (int i = 0; i < fe.dofs_per_cell; ++i) {
    for (int j = 0; j < fe.dofs_per_cell; ++j) {
      cell_matrix(i,j) = 1;
    }
  }
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

  dealii::PETScWrappers::MPI::SparseMatrix index_hits;
  index_hits.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  for (auto cell : cells) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      for (auto index_j : local_dof_indices) {
        index_hits.add(index_i, index_j, 1);
      }
    }
    index_hits.compress(dealii::VectorOperation::add);
    system_matrix.add(local_dof_indices, local_dof_indices, cell_matrix);
    system_matrix.compress(dealii::VectorOperation::add);
  }

  for (int i = 0; i < local_dof_indices.size(); ++i) {
    for (int j = 0; j < local_dof_indices.size(); ++j) {
      EXPECT_EQ(system_matrix(local_dof_indices[i], local_dof_indices[j]),
                index_hits(local_dof_indices[i], local_dof_indices[j]));
    }
  }
}
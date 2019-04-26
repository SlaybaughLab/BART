#ifndef BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_
#define BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_

#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include <gtest/gtest.h>

namespace bart {

namespace testing {

template <int dim>
struct TriangulationType {
  using type = dealii::parallel::distributed::Triangulation<dim>;
};

template <>
struct TriangulationType<1> {
  using type = dealii::Triangulation<1>;
};

template struct TriangulationType<2>;
template struct TriangulationType<3>;

template <int dim>
 class DealiiTestDomain {
  public:
   static constexpr int dimension = dim;
   DealiiTestDomain();
   void SetUpDealii();

   using Cell = typename dealii::DoFHandler<dim>::active_cell_iterator;
   dealii::ConstraintMatrix constraint_matrix_;
   typename TriangulationType<dim>::type triangulation_;
   dealii::DoFHandler<dim> dof_handler_;
   dealii::FE_Q<dim> fe_;
   dealii::IndexSet locally_relevant_dofs;
   dealii::IndexSet locally_owned_dofs_;
   std::vector<Cell> cells_;
   dealii::DynamicSparsityPattern dsp_;

   dealii::PETScWrappers::MPI::SparseMatrix matrix_1, matrix_2, matrix_3;

  private:
   void SetUpDofs();
};

template <int dim>
inline DealiiTestDomain<dim>::DealiiTestDomain()
    : triangulation_(MPI_COMM_WORLD,
                     typename dealii::Triangulation<dim>::MeshSmoothing(
                         dealii::Triangulation<dim>::smoothing_on_refinement |
                             dealii::Triangulation<dim>::smoothing_on_coarsening)),
      dof_handler_(triangulation_),
      fe_(1) {}

template <>
inline DealiiTestDomain<1>::DealiiTestDomain()
    : triangulation_(typename dealii::Triangulation<1>::MeshSmoothing(
    dealii::Triangulation<1>::smoothing_on_refinement |
        dealii::Triangulation<1>::smoothing_on_coarsening)),
      dof_handler_(triangulation_),
      fe_(1) {}

template <int dim>
inline void DealiiTestDomain<dim>::SetUpDealii() {
  dealii::GridGenerator::hyper_cube(triangulation_, 0, 1);

  if (dim == 1)
    triangulation_.refine_global(4);
  else
    triangulation_.refine_global(2);

  SetUpDofs();

  for (auto cell = dof_handler_.begin_active(); cell != dof_handler_.end(); ++ cell) {
    if (cell->is_locally_owned())
      cells_.push_back(cell);
  }

  matrix_1.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp_, MPI_COMM_WORLD);
  matrix_2.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp_, MPI_COMM_WORLD);
  matrix_3.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp_, MPI_COMM_WORLD);
}

template <int dim>
inline void DealiiTestDomain<dim>::SetUpDofs() {
  dof_handler_.distribute_dofs(fe_);
  locally_owned_dofs_ = dof_handler_.locally_owned_dofs();
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_,
                                                  locally_relevant_dofs);

  constraint_matrix_.clear();
  constraint_matrix_.reinit(locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();

  dsp_.reinit(locally_relevant_dofs.size(),
              locally_relevant_dofs.size(),
              locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp_,
                                          constraint_matrix_, false);

  dealii::SparsityTools::distribute_sparsity_pattern(
      dsp_,
      dof_handler_.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs);

  constraint_matrix_.condense(dsp_);
}

template <>
inline void DealiiTestDomain<1>::SetUpDofs() {
  auto n_mpi_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  dealii::GridTools::partition_triangulation(n_mpi_processes, triangulation_);
  dof_handler_.distribute_dofs(fe_);
  dealii::DoFRenumbering::subdomain_wise(dof_handler_);

  auto locally_owned_dofs_vector =
      dealii::DoFTools::locally_owned_dofs_per_subdomain(dof_handler_);

  locally_owned_dofs_ = locally_owned_dofs_vector.at(this_process);

  constraint_matrix_.clear();
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();

  dsp_.reinit(dof_handler_.n_dofs(), dof_handler_.n_dofs());
  dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp_,
                                          constraint_matrix_, false);
}

using DealiiTestDomains = ::testing::Types<bart::testing::DealiiTestDomain<1>,
                                            bart::testing::DealiiTestDomain<2>,
                                            bart::testing::DealiiTestDomain<3>>;

} // namespace testing

} // namespace bart

#endif // BART_SRC_TEST_HELPERS_MPI_TEST_FIXTURE_H_

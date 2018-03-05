#include "../../src/mesh/mesh_generator.h"
#include "../test_utilities.h"

#include <deal.II/base/types.h>

template <int dim>
void SetupParameters (dealii::ParameterHandler &prm)
{
  prm.declare_entry ("is mesh generated by deal.II", "true",
                     dealii::Patterns::Bool(), "");
  prm.declare_entry ("have reflective BC", "true",
                     dealii::Patterns::Bool(), "");
  prm.declare_entry ("uniform refinements", "1",
                     dealii::Patterns::Integer(), "");
  prm.declare_entry ("x, y, z max values of boundary locations", "1.0,1.0,1.0",
                     dealii::Patterns::List (dealii::Patterns::Double ()), "");
  prm.declare_entry ("number of cells for x, y, z directions", "2,2,2",
                     dealii::Patterns::List (dealii::Patterns::Integer ()), "");
  prm.declare_entry ("reflective boundary names", "xmin, ymax",
                     dealii::Patterns::List(dealii::Patterns::Anything ()), "");
  prm.enter_subsection ("material ID map");
  {
    std::string id_fname = SOURCE_DIR + std::string ("/matid.homogeneous.")
        + std::to_string (dim) + std::string ("d");
    prm.declare_entry ("material id file name", id_fname,
                       dealii::Patterns::FileName(), "file name for material id map");
  }
  prm.leave_subsection ();
}

template <int dim>
void Test (dealii::ParameterHandler &prm)
{
  dealii::deallog.push (dealii::Utilities::int_to_string(dim)+"D");
  
  auto process_id = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  auto n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  
  AssertThrow (!(n_proc&(n_proc-1)),
               dealii::ExcMessage("Number of cores should be power of 2."));
  
  // triangulation for the grid in the scope of the test function
  dealii::parallel::distributed::Triangulation<dim> tria (MPI_COMM_WORLD);
  
  // make grid
  std::unique_ptr<MeshGenerator<dim>> msh_ptr =
      std::unique_ptr<MeshGenerator<dim>> (new MeshGenerator<dim>(prm));
  msh_ptr->MakeGrid (tria);
  
  AssertThrow (tria.n_locally_owned_active_cells()==std::pow(2, 2*dim)/n_proc,
               dealii::ExcInternalError());
  
  tria.refine_global (2);
  
  AssertThrow (tria.n_locally_owned_active_cells()==std::pow(2, 4*dim)/n_proc,
               dealii::ExcInternalError());
  
  dealii::deallog << "Global refinements check OK." << std::endl;
  
  for (typename dealii::Triangulation<dim>::active_cell_iterator
       cell=tria.begin_active(); cell!=tria.end(); ++cell)
    if (cell->is_locally_owned())
      // material ID should be input - 1 = 111 - 1 = 110
      AssertThrow (cell->material_id()==110, dealii::ExcInternalError());
  
  dealii::deallog << "Material ID check OK." << std::endl;
  
  dealii::deallog.pop ();
}

int main (int argc, char *argv[])
{
  // initialize MPI and log and declare ParameterHandler object
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  testing::MPILogInit init;
  dealii::ParameterHandler prm;
  
  // parameter processing
  SetupParameters<2> (prm);
  
  // testing 2D
  Test<2> (prm);
  
  // clearing prm so new parameters can be set for other dimensions if needed
  prm.clear ();
  testing::deallogstream << std::endl;
  
  // 3D testing section
  SetupParameters<3> (prm);
  Test<3> (prm);
  prm.clear ();
  
  return 0;
}

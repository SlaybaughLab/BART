#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "equation_base.h"
#include "../aqdata/aq_base.h"
#include "../aqdata/aq_lsgc.h"
#include "../common/bart_tools.h"

using namespace dealii;

template <int dim>
BartDriver<dim>::BartDriver (ParameterHandler &prm)
:
triangulation (MPI_COMM_WORLD,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
transport_model_name(prm.get("transport model")),
aq_name(prm.get("angular quadrature name")),
n_group(prm.get_integer("number of groups")),
n_azi(prm.get_integer("angular quadrature order")),
is_eigen_problem(prm.get_bool("do eigenvalue calculations")),
do_nda(prm.get_bool("do NDA")),
do_print_sn_quad(prm.get_bool("do print angular quadrature info")),
have_reflective_bc(prm.get_bool("have reflective BC")),
p_order(prm.get_integer("finite element polynomial degree")),
global_refinements(prm.get_integer("uniform refinements")),
namebase(prm.get("output file name base")),
ho_linear_solver_name(prm.get("HO linear solver name")),
ho_preconditioner_name(prm.get("HO preconditioner name"))
{
  
  aqd_ptr = build_aq_model (prm)
  n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
  n_azi = aqd_ptr->get_sn_order ();
  n_dir = aqd_ptr->get_n_dir ();
  msh_ptr = build_mesh (prm);
  mat_ptr = build_material (prm);
  fe = build_finite_element (prm);
  sflxes_proc.resize (n_group);
}

template <int dim>
BartDriver<dim>::~BartDriver ()
{
  dof_handler.clear();
}

template <int dim>
void BartDriver<dim>::report_system ()
{
  pout << "SN quadrature order: " << n_azi << std::endl
  << "Number of angles: " << n_dir << std::endl
  << "Number of groups: " << n_group << std::endl;

  pout << "Transport model: " << transport_model_name << std::endl;
  pout << "Spatial discretization: " << discretization << std::endl;
  pout << "HO linear solver: " << ho_linear_solver_name << std::endl;
  if (ho_linear_solver_name!="direct")
    pout << "HO preconditioner: " << ho_preconditioner_name << std::endl;
  pout << "do NDA? " << do_nda << std::endl;
  
  pout << "Number of cells: " << triangulation.n_global_active_cells() << std::endl;
  pout << "High-order total DoF counts: " << n_total_ho_vars*dof_handler.n_dofs() << std::endl;

  if (is_eigen_problem)
    pout << "Problem type: k-eigenvalue problem" << std::endl;
  if (do_nda)
    pout << "NDA total DoF counts: " << n_group*dof_handler.n_dofs()) << std::endl;
  pout << "print sn quad? " << do_print_sn_quad << std::endl;
  pout << "is eigenvalue problem? " << is_eigen_problem << std::endl;
}

template <int dim>
void BartDriver<dim>::setup_system ()
{
  radio ("setup system");
  dof_handler.distribute_dofs (*fe);
  
  local_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           relevant_dofs);
  
  constraints.clear ();
  constraints.reinit (relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();
  
  DynamicSparsityPattern dsp (relevant_dofs);
  
  if (discretization=="dfem")
    DoFTools::make_flux_sparsity_pattern (dof_handler,
                                          dsp,
                                          constraints,
                                          false);
  else
    DoFTools::make_sparsity_pattern (dof_handler,
                                     dsp,
                                     constraints,
                                     false);
  
  // setting up dsp with telling communicator and relevant dofs
  SparsityTools::distribute_sparsity_pattern
  (dsp,
   dof_handler.n_locally_owned_dofs_per_processor (),
   MPI_COMM_WORLD,
   relevant_dofs);
  
  itr_ptr = build_iterations (prm, dof_handler, msh_ptr, aqd_ptr, mat_ptr);
  itr_ptr->initialize_system_matrices_vectors (dsp, local_dofs);
}

template <int dim>
void BartDriver<dim>::output_results () const
{
  std::string sec_name = "Graphical output";
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  if (is_eigen_problem)
    data_out.add_data_vector (keff, "keff");
  
  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "phi_g_" << g;
    data_out.add_data_vector (sflxes_proc[g], os.str ());
  }

  Vector<float> subdomain (triangulation.n_active_cells ());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain ();
  data_out.add_data_vector (subdomain, "subdomain");
  data_out.build_patches ();

  const std::string filename =
  (namebase + "-" + discretization + "-" + Utilities::int_to_string
   (triangulation.locally_owned_subdomain (), 4));
  std::ofstream output ((filename + ".vtu").c_str ());
  data_out.write_vtu (output);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
         ++i)
      filenames.push_back (namebase + "-" + discretization + "-" +
                           Utilities::int_to_string (i) + ".vtu");
    std::ostringstream os;
    os << namebase << "-" << discretization << "-" << global_refinements << ".pvtu";
    std::ofstream master_output ((os.str()).c_str ());
    data_out.write_pvtu_record (master_output, filenames);
  }
}

template <int dim>
void BartDriver<dim>::run ()
{
  msh_ptr->make_grid (triangulation);
  setup_system ();
  report_system ();
  // solve the problem using iterative methods specified in Iterations class
  itr_ptr->solve_problems (sflxes_proc);
  if (is_eigen_problem)
    itr_ptr->get_keff (keff);
  output_results ();
}

// explicit instantiation to avoid linking error
template class BartDriver<2>;
template class BartDriver<3>;

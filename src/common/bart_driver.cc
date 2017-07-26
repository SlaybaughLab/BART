#include <deal.II/fe/fe_values.h>

#include <boost/algorithm/string.hpp>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/cell_id.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <algorithm>

#include "transport_base.h"
#include "../aqdata/aq_base.h"
#include "../aqdata/aq_lsgc.h"

using namespace dealii;

template <int dim>
BartDriver<dim>::BartDriver (ParameterHandler &prm)
:
mpi_communicator (MPI_COMM_WORLD),
triangulation (mpi_communicator,
               typename Triangulation<dim>::MeshSmoothing
               (Triangulation<dim>::smoothing_on_refinement |
                Triangulation<dim>::smoothing_on_coarsening)),
dof_handler (triangulation),
prm(prm),
transport_model_name(prm.get("transport model")),
aq_name(prm.get("angular quadrature name")),
discretization(prm.get("spatial discretization")),
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
ho_preconditioner_name(prm.get("HO preconditioner name")),
pcout(std::cout,
      (Utilities::MPI::this_mpi_process(mpi_communicator)
       == 0))
{
  aqd_ptr = build_aq_model (prm)
  aqd_ptr->make_aq (prm);
  n_total_ho_vars = aqd_ptr->get_n_total_ho_vars ();
  n_azi = aqd_ptr->get_sn_order ();
  n_dir = aqd_ptr->get_n_dir ();
  msh_ptr = std_cxx11::shared_ptr<MeshGenerator<dim> >
  (new MeshGenerator<dim>(prm));
  this->process_input ();
}

template <int dim>
BartDriver<dim>::~BartDriver ()
{
  dof_handler.clear();
}

template <int dim>
void BartDriver<dim>::report_system ()
{
  pcout << "SN quadrature order: " << n_azi << std::endl
  << "Number of angles: " << n_dir << std::endl
  << "Number of groups: " << n_group << std::endl;

  radio ("Transport model", transport_model_name);
  radio ("Spatial discretization", discretization);
  radio ("HO linear solver", ho_linear_solver_name);
  if (ho_linear_solver_name!="direct")
    radio ("HO preconditioner", ho_preconditioner_name);
  radio ("do NDA?", do_nda);
  
  radio ("Number of cells", triangulation.n_global_active_cells());
  radio ("High-order total DoF counts", n_total_ho_vars*dof_handler.n_dofs());

  if (is_eigen_problem)
    radio ("Problem type: k-eigenvalue problem");
  if (do_nda)
    radio ("NDA total DoF counts", n_group*dof_handler.n_dofs());
  radio ("print sn quad?", do_print_sn_quad);
  if (do_print_sn_quad &&
      Utilities::MPI::this_mpi_process(mpi_communicator)==0)
    aqd_ptr->print_angular_quad ();
  radio ("is eigenvalue problem?", is_eigen_problem);
}

template <int dim>
void BartDriver<dim>::setup_system ()
{
  radio ("setup system");
  initialize_dealii_objects ();
  initialize_system_matrices_vectors ();
}

template <int dim>
void BartDriver<dim>::initialize_system_matrices_vectors ()
{
  DynamicSparsityPattern dsp (relevant_dofs);

  if (discretization=="dfem")
  {
    /*
     Table<2,DoFTools::Coupling> cell_coupling (1,1);
     Table<2,DoFTools::Coupling> face_coupling (1,1);

     cell_coupling[0][0] = DoFTools::nonzero;
     face_coupling[0][0] = DoFTools::nonzero;

     DoFTools::make_flux_sparsity_pattern (dof_handler,
     dsp,
     cell_coupling,
     face_coupling);
     */

    DoFTools::make_flux_sparsity_pattern (dof_handler,
                                          dsp,
                                          constraints,
                                          false);
  }
  else
    DoFTools::make_sparsity_pattern (dof_handler,
                                     dsp,
                                     constraints,
                                     false);

  // be careful with the following
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.n_locally_owned_dofs_per_processor (),
                                              mpi_communicator,
                                              relevant_dofs);

  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys.push_back (new LA::MPI::SparseMatrix);
      vec_lo_rhs.push_back (new LA::MPI::Vector);
      vec_lo_sflx.push_back (new LA::MPI::Vector);
      vec_lo_sflx_old.push_back (new LA::MPI::Vector);
      vec_lo_fixed_rhs.push_back (new LA::MPI::Vector);
    }

    vec_ho_sflx.push_back (new LA::MPI::Vector);
    vec_ho_sflx_prev_gen.push_back (new LA::MPI::Vector);
    vec_ho_sflx_old.push_back (new LA::MPI::Vector);

    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys.push_back (new LA::MPI::SparseMatrix);
      vec_aflx.push_back (new LA::MPI::Vector);
      vec_ho_rhs.push_back (new LA::MPI::Vector);
      vec_ho_fixed_rhs.push_back (new LA::MPI::Vector);
    }
  }

  for (unsigned int g=0; g<n_group; ++g)
  {
    if (do_nda)
    {
      vec_lo_sys[g]->reinit (local_dofs,
                             local_dofs,
                             dsp,
                             mpi_communicator);
      vec_lo_rhs[g]->reinit (local_dofs,
                             mpi_communicator);
      vec_lo_fixed_rhs[g]->reinit (local_dofs,
                                   mpi_communicator);
      vec_lo_sflx[g]->reinit (local_dofs,
                              mpi_communicator);
      vec_lo_sflx_old[g]->reinit (local_dofs,
                                  mpi_communicator);
    }

    vec_ho_sflx[g]->reinit (local_dofs,
                            mpi_communicator);
    vec_ho_sflx_old[g]->reinit (local_dofs,
                                mpi_communicator);

    for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
    {
      vec_ho_sys[get_component_index(i_dir, g)]->reinit(local_dofs,
                                                        local_dofs,
                                                        dsp,
                                                        mpi_communicator);
      vec_aflx[get_component_index(i_dir, g)]->reinit(local_dofs,
                                                      mpi_communicator);
      vec_ho_rhs[get_component_index(i_dir, g)]->reinit (local_dofs,
                                                         mpi_communicator);
      vec_ho_fixed_rhs[get_component_index(i_dir, g)]->reinit (local_dofs,
                                                               mpi_communicator);
    }
  }
}

template <int dim>
void BartDriver<dim>::initialize_dealii_objects ()
{
  if (discretization=="dfem")
    fe = (new FE_DGQ<dim> (p_order));
  else
    fe = (new FE_Q<dim> (p_order));

  dof_handler.distribute_dofs (*fe);

  local_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler,
                                           relevant_dofs);

  constraints.clear ();
  constraints.reinit (relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  constraints.close ();
}

template <int dim>
void BartDriver<dim>::output_results () const
{
  std::string sec_name = "Graphical output";
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);

  for (unsigned int g=0; g<n_group; ++g)
  {
    std::ostringstream os;
    os << "ho_phi_g_" << g;
    data_out.add_data_vector (sflx_proc[g], os.str ());
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

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i=0;
         i<Utilities::MPI::n_mpi_processes(mpi_communicator);
         ++i)
      filenames.push_back (namebase + "-" + discretization + "-" +
                           Utilities::int_to_string (i, 4) + ".vtu");
    std::ostringstream os;
    os << namebase << "-" << discretization << "-" << global_refinements << ".pvtu";
    std::ofstream master_output ((os.str()).c_str ());
    data_out.write_pvtu_record (master_output, filenames);
  }
}

template <int dim>
void BartDriver<dim>::run ()
{
  radio ("making grid");
  msh_ptr->make_grid (triangulation);
  msh_ptr->get_relevant_cell_iterators (dof_handler,
                                        local_cells,
                                        ref_bd_cells,
                                        is_cell_at_bd,
                                        is_cell_at_ref_bd);
  setup_system ();
  report_system ();
  assemble_ho_system ();
  do_iterations ();
  output_results();
}

// wrapper functions used to retrieve info from various Hash tables
template <int dim>
unsigned int BartDriver<dim>::get_component_index
(unsigned int incident_angle_index, unsigned int g)
{
  // retrieve component indecis given direction and group
  // must be used after initializing the index map
  return component_index[std::make_pair (incident_angle_index, g)];
}

template <int dim>
unsigned int BartDriver<dim>::get_component_direction (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].first;
}

template <int dim>
unsigned int BartDriver<dim>::get_component_group (unsigned int comp_ind)
{
  return inverse_component_index[comp_ind].second;
}

template <int dim>
unsigned int BartDriver<dim>::get_reflective_direction_index
(unsigned int boundary_id, unsigned int incident_angle_index)
{
  AssertThrow (is_reflective_bc[boundary_id],
               ExcMessage ("must be reflective boundary to retrieve the reflective boundary"));
  return reflective_direction_index[std::make_pair (boundary_id,
                                                    incident_angle_index)];
}

//functions used to cout information for diagonose or just simply cout
template <int dim>
void BartDriver<dim>::radio (std::string str)
{
  pcout << str << std::endl;
}

template <int dim>
void BartDriver<dim>::radio (std::string str1, std::string str2)
{
  pcout << str1 << ": " << str2 << std::endl;
}

template <int dim>
void BartDriver<dim>::radio (std::string str,
                                double num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void BartDriver<dim>::radio (std::string str1, unsigned int num1,
                                std::string str2, unsigned int num2,
                                std::string str3, unsigned int num3)
{
  pcout << str1 << ": " << num1 << ", ";
  pcout << str2 << ": " << num2 << ", ";
  pcout << str3 << ": " << num3 << std::endl;;
}

template <int dim>
void BartDriver<dim>::radio (std::string str,
                                unsigned int num)
{
  pcout << str << ": " << num << std::endl;
}

template <int dim>
void BartDriver<dim>::radio (std::string str, bool boolean)
{
  pcout << str << ": " << (boolean?"true":"false") << std::endl;
}

template <int dim>
void BartDriver<dim>::radio ()
{
  pcout << "-------------------------------------" << std::endl << std::endl;
}

// explicit instantiation to avoid linking error
template class BartDriver<2>;
template class BartDriver<3>;

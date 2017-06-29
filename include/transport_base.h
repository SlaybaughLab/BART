#ifndef __transport_base_h__
#define __transport_base_h__
#include <deal.II/lac/generic_linear_algebra.h>
namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
}

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_poly.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "problem_definition.h"
#include "mesh_generator.h"
#include "material_properties.h"
#include "aq_base.h"

using namespace dealii;

template <int dim>
class TransportBase
{
public:
  TransportBase (ParameterHandler &prm);// : ProblemDefinition<dim> (prm){}
  virtual ~TransportBase ();
  
  void run ();
  
  virtual void pre_assemble_cell_matrices
  (const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  virtual void integrate_cell_bilinear_form
  (const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   unsigned int &i_dir,
   unsigned int &g,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  virtual void integrate_boundary_bilinear_form
  (const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   unsigned int &i_dir,
   unsigned int &g);
  
  virtual void integrate_interface_bilinear_form
  (const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
   unsigned int &fn,/*concerning face number in local cell*/
   unsigned int &i_dir,
   unsigned int &g,
   FullMatrix<double> &vp_up,
   FullMatrix<double> &vp_un,
   FullMatrix<double> &vn_up,
   FullMatrix<double> &vn_un);
  
  virtual void generate_moments ();
  virtual void postprocess ();
  virtual void generate_ho_source ();
  
private:
  void setup_system ();
  void generate_globally_refined_grid ();
  void report_system ();
  void print_angular_quad ();
  
  // void setup_lo_system();
  void setup_boundary_ids ();
  void get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                      std::vector<double> &local_mfps);
  void assemble_ho_volume_boundary ();
  void assemble_ho_interface ();
  void assemble_ho_system ();
  void do_iterations ();
  void process_input ();
  void initialize_material_id ();
  void initialize_dealii_objects ();
  void initialize_system_matrices_vectors ();
  unsigned int get_component_index (unsigned int &incident_angle_index, unsigned int &g);
  unsigned int get_direction (unsigned int &comp_ind);
  unsigned int get_component_group (unsigned int &comp_ind);
  
  void assemble_lo_system ();
  void prepare_correction_aflx ();
  
  void initialize_ho_preconditioners ();
  void ho_solve ();
  void lo_solve ();
  void refine_grid ();
  void output_results () const;
  void generate_fixed_source ();
  void power_iteration ();
  void source_iteration ();
  void scale_fiss_transfer_matrices ();
  void renormalize_sflx (std::vector<LA::MPI::Vector*> &target_sflxes,
                         double normalization_factor);
  
  void radio (std::string str);
  void radio (std::string str1, std::string str2);
  void radio (std::string str, double num);
  void radio (std::string str, unsigned int num);
  
  double estimate_k (double &fiss_source,
                     double &fiss_source_prev_gen,
                     double &k_prev_gen);
  double estimate_fiss_source (std::vector<Vector<double> > &phis_this_process);
  double estimate_phi_diff (std::vector<LA::MPI::Vector*> &phis_newer,
                            std::vector<LA::MPI::Vector*> &phis_older);
  //Properties<dim> *mat_prop;
  
  void NDA_PI ();
  void NDA_SI ();
  void initialize_aq (ParameterHandler &prm);
  
  std_cxx11::shared_ptr<ProblemDefinition> def_ptr;
  std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr;
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;
  std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr;
  
  std_cxx11::shared_ptr<FEValues<dim> > fv;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei;
  
  std::string transport_model_name;
  
  std::vector<typename DoFHandler<dim>::active_cell_iterator> local_cells;
  std::vector<bool> is_cell_at_bd;
  std::vector<bool> is_cell_at_ref_bd;
protected:
  unsigned int get_reflective_direction_index (unsigned int &boundary_id,
                                               unsigned int &incident_angle_index);
  
  //std_cxx11::shared_ptr<FE_Poly<TensorProductPolynomials<dim>,dim,dim> > fe;
  
  FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe;
  std_cxx11::shared_ptr<QGauss<dim> > q_rule;
  std_cxx11::shared_ptr<QGauss<dim-1> > qf_rule;
  
  
  MPI_Comm mpi_communicator;
  
  parallel::distributed::Triangulation<dim> triangulation;
  
  DoFHandler<dim> dof_handler;
  // FE_DGQ<dim> *fe;
  
  // FixIt: involve relevant_dofs for future if refinement is necessary
  IndexSet local_dofs;
  IndexSet relevant_dofs;
  
  const double err_k_tol;
  const double err_phi_tol;
  
  double k_ho;
  double k_ho_prev_gen;
  double total_angle;
  double c_penalty;
  double fission_source;
  double fission_source_prev_gen;
  
  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  bool is_explicit_reflective;
  bool do_print_sn_quad;
  
  unsigned int n_q;
  unsigned int n_qf;
  unsigned int dofs_per_cell;
  
  unsigned int n_dir;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  unsigned int n_group;
  unsigned int n_material;
  unsigned int p_order;
  unsigned int global_refinements;
  
  std::string discretization;
  std::string namebase;
  std::string aq_name;
  
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> neigh_dof_indices;
  
  // HO system
  std::vector<LA::MPI::SparseMatrix*> vec_ho_sys;
  std::vector<LA::MPI::Vector*> vec_aflx;
  std::vector<LA::MPI::Vector*> vec_ho_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_ho_sflx;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_old;
  std::vector<LA::MPI::Vector*> vec_ho_sflx_prev_gen;
  
  // LO system
  std::vector<LA::MPI::SparseMatrix*> vec_lo_sys;
  std::vector<LA::MPI::Vector*> vec_lo_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_fixed_rhs;
  std::vector<LA::MPI::Vector*> vec_lo_sflx;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_old;
  std::vector<LA::MPI::Vector*> vec_lo_sflx_prev_gen;
  
  std::vector<Tensor<1, dim> > omega_i;
  std::vector<double> wi;
  std::vector<double> tensor_norms;
  std::vector<std::vector<double> > all_sigt;
  std::vector<std::vector<double> > all_inv_sigt;
  std::vector<std::vector<double> > all_q;
  std::vector<std::vector<double> > all_q_per_ster;
  std::vector<std::vector<double> > all_nusigf;
  std::vector<std::vector<std::vector<double> > > all_sigs;
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf_per_ster;
  std::vector<std::vector<std::vector<double> > > ho_scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > lo_scaled_fiss_transfer;
  std::vector<Vector<double> > sflx_this_processor;
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::set<unsigned int> fissile_ids;
  
  ConditionalOStream pcout;
  
  std::vector<std_cxx11::shared_ptr<LA::MPI::PreconditionAMG> > pre_ho_amg;
  
  ConstraintMatrix constraints;
};

#endif	// define  __transport_base_h__

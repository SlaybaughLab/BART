#ifndef __transport_base_h__
#define __transport_base_h__
#include <deal.II/lac/generic_linear_algebra.h>

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

#include "../common/problem_definition.h"
#include "../common/preconditioner_solver.h"
#include "../mesh/mesh_generator.h"
#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"

using namespace dealii;

template <int dim>
class EquationBase
{
public:
  // For SN and NDA (NDA still needs quadrature to calculate corrections)
  EquationBase (ParameterHandler &prm,
                 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
  
  // For diffusion (future work for someone)
  EquationBase (ParameterHandler &prm,
                const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
  
  virtual ~EquationBase ();
  
  void run ();
  
  virtual void assemble_volume_boundary
  (std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
   std::vector<bool> &is_cell_at_bd);
  virtual void assemble_interface
  (std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells);
  virtual void assemble_system
  (std::vector<typename DoFHandler<dim>::active_cell_iterator> &local_cells,
   std::vector<bool> &is_cell_at_bd);
  
  virtual void pre_assemble_cell_matrices
  (const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp);
  
  virtual void integrate_cell_bilinear_form
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   const std_cxx11::shared_ptr<FEValues<dim> > fv,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   FullMatrix<double> &cell_matrix,
   std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
   std::vector<FullMatrix<double> > &collision_at_qp,
   const unsigned int &g,
   const unsigned int &i_dir=0);
  
  virtual void integrate_boundary_bilinear_form
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   FullMatrix<double> &cell_matrix,
   const unsigned int &g,
   const unsigned int &i_dir=0);
  
  virtual void integrate_reflective_boundary_linear_form
  (const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
   typename DoFHandler<dim>::active_cell_iterator &cell,
   unsigned int &fn,/*face number*/
   std::vector<Vector<double> > &cell_rhses,
   const unsigned int &g,
   const unsigned int &i_dir=0);
  
  virtual void integrate_interface_bilinear_form
  (std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
   const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
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
  
  virtual void generate_moments
  (std::vector<Vector<double>*> &vec_ho_sflx,
   std::vector<Vector<double>*> &vec_ho_sflx_old,
   std::vector<Vector<double>*> &sflx_proc);
  virtual void postprocess ();
  virtual void generate_ho_rhs ();
  virtual void generate_ho_fixed_source
  (std::vector<PETScWrappers::MPI::Vector*> &vec_ho_fixed_rhs,
   std::vector<Vector<double> > &sflx_this_proc);
  
  // override this three functions in derived classes
  virtual unsigned int get_component_index (unsigned int incident_angle_index, unsigned int g);
  virtual unsigned int get_component_direction (unsigned int comp_ind);
  virtual unsigned int get_component_group (unsigned int comp_ind);
  
private:
  void setup_system ();
  void generate_globally_refined_grid ();
  void report_system ();
  void print_angular_quad ();
  
  // void setup_lo_system();
  void setup_boundary_ids ();
  void get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                      std::vector<double> &local_mfps);
  
  void process_input (const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                      const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                      const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);
  void initialize_system_matrices_vectors ();
  void prepare_correction_aflx ();
  void initialize_fiss_process ();
  void update_ho_moments_in_fiss ();
  void update_fiss_source_keff ();
  void source_iteration ();
  void scale_fiss_transfer_matrices ();
  void renormalize_sflx (std::vector<LA::MPI::Vector*> &target_sflxes);
  void initialize_aq (ParameterHandler &prm);
  void initialize_assembly_related_objects
  (FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe);
  
  double estimate_fiss_source (std::vector<Vector<double> > &phis_this_process);
  
  std_cxx11::shared_ptr<MaterialProperties> mat_ptr;
  std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr;
  std_cxx11::shared_ptr<PreconditionerSolver> sol_ptr;
  
  std::string transport_model_name;
  std::string ho_linear_solver_name;
  std::string ho_preconditioner_name;
  std::string discretization;
  std::string namebase;
  std::string aq_name;
  
protected:
  unsigned int get_reflective_direction_index (unsigned int boundary_id,
                                               unsigned int incident_angle_index);
  
  // "c" in the following quantities means "correction" for NDA use
  std_cxx11::shared_ptr<QGauss<dim> > q_rule;
  std_cxx11::shared_ptr<QGauss<dim-1> > qf_rule;
  std_cxx11::shared_ptr<QGauss<dim> > qc_rule;
  std_cxx11::shared_ptr<QGauss<dim-1> > qfc_rule;
  std_cxx11::shared_ptr<FEValues<dim> > fv;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei;
  std_cxx11::shared_ptr<FEValues<dim> > fvc;
  std_cxx11::shared_ptr<FEFaceValues<dim> > fvfc;
  
  double keff;
  double keff_prev_gen;
  double fission_source;
  double fission_source_prev_gen;
  
  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  
  const unsigned int nda_quadrature_order;
  unsigned int n_q;
  unsigned int n_qf;
  unsigned int n_qc;
  unsigned int n_qfc;
  unsigned int dofs_per_cell;
  
  unsigned int n_dir;
  unsigned int n_azi;
  unsigned int n_total_vars;
  unsigned int n_group;
  unsigned int n_material;
  unsigned int p_order;
  unsigned int global_refinements;
  
  std::vector<unsigned int> linear_iters;
  
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<types::global_dof_index> neigh_dof_indices;
  
  
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
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > scat_scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > scaled_fiss_transfer;
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::set<unsigned int> fissile_ids;
};

#endif	// define  __transport_base_h__

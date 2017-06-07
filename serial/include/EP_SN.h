#ifndef __ep_sn_h__
#define __ep_sn_h__
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/sparse_matrix.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_poly.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "problem_definition.h"

using namespace dealii;

template <int dim>
class EP_SN : public ProblemDefinition<dim>
{
public:
  EP_SN (ParameterHandler &prm);// : ProblemDefinition<dim> (prm){}
  ~EP_SN ();

  void run ();

private:
  void setup_system ();
  void generate_globally_refined_grid ();
  void report_system ();
  // void setup_lo_system();
  void setup_boundary_ids ();
  void get_cell_mfps (unsigned int &material_id, double &cell_dimension,
                      std::vector<double> &local_mfps);
  void assemble_ho_volume_boundary ();
  void assemble_ho_volume_boundary_full ();
  void assemble_ho_volume_boundary_third ();
  void assemble_ho_interface ();
  void assemble_ho_system ();
  void do_iterations ();
  void process_input ();
  void initialize_material_id ();
  unsigned int get_component_index (unsigned int &incident_angle_index, unsigned int &g);
  unsigned int get_direction (unsigned int &comp_ind);
  unsigned int get_component_group (unsigned int &comp_ind);
  unsigned int get_reflective_direction_index (unsigned int &boundary_id,
                                               unsigned int &incident_angle_index);
  void get_cell_relative_position (Point<dim> &position,
                                   std::vector<unsigned int> &relative_position);

  void assemble_lo_system ();
  void prepare_correction_aflx ();

  void initialize_ho_preconditioners ();
  void ho_solve ();
  void lo_solve ();
  void refine_grid ();
  void output_results () const;
  void generate_moments ();
  void generate_ho_source ();
  void generate_ho_source_new ();
  void generate_fixed_source ();
  void generate_fixed_source_new ();
  void power_iteration ();
  void source_iteration ();
  void scale_fiss_transfer_matrices ();
  void renormalize_sflx (std::vector<Vector<double>*> &target_sflxes, double &normalization_factor);

  void local_matrix_check (FullMatrix<double> &local_mat,
                           std::string str,
                           unsigned int ind);

  void local_radio (std::string str);

  void global_matrix_check (unsigned int ind);

  double estimate_k (double &fiss_source,
                     double &fiss_source_prev_gen,
                     double &k_prev_gen);
  double estimate_fiss_source (std::vector<Vector<double>*> &phis);
  double estimate_phi_diff (std::vector<Vector<double>*> &phis_newer,
                            std::vector<Vector<double>*> &phis_older);
  //Properties<dim> *mat_prop;

  void NDA_PI ();
  void NDA_SI ();

  ProblemDefinition<dim>* paras;

  Triangulation<dim> triangulation;

  DoFHandler<dim> dof_handler;
  // FE_DGQ<dim> *fe;
  FE_Poly<TensorProductPolynomials<dim>,dim,dim> *fe;
  
  double pi;

  // HO system
  std::vector<SparseMatrix<double>*> vec_ho_sys;
  std::vector<Vector<double>*> vec_aflx;
  std::vector<Vector<double>*> vec_ho_rhs;
  std::vector<Vector<double>*> vec_ho_fixed_rhs;
  std::vector<Vector<double>*> vec_ho_sflx;
  std::vector<Vector<double>*> vec_ho_sflx_old;
  std::vector<Vector<double>*> vec_ho_sflx_prev_gen;
  double k_ho;
  double k_ho_prev_gen;

  // LO system
  std::vector<SparseMatrix<double>*> vec_lo_sys;
  std::vector<Vector<double>*> vec_lo_rhs;
  std::vector<Vector<double>*> vec_lo_fixed_rhs;
  std::vector<Vector<double>*> vec_lo_sflx;
  std::vector<Vector<double>*> vec_lo_sflx_old;
  std::vector<Vector<double>*> vec_lo_sflx_prev_gen;

  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;

  std::string discretization;
  double total_angle;
  unsigned int n_azi;
  unsigned int n_total_ho_vars;
  unsigned int n_group;
  unsigned int n_material;
  int p_order;
  int global_refinements;
  double c_penalty;
  unsigned int n_dir;
  std::set<unsigned int> fissile_ids;
  std::vector<unsigned int> ncell_per_dir;
  std::vector<double> cell_size_all_dir;
  std::vector<double> axis_max_values;

  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  bool is_explicit_reflective;
  bool do_print_sn_quad;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  std::unordered_map<unsigned int, bool> is_material_fissile;

  std::vector<Tensor<1, dim> > omega_i;
  std::vector<double> wi;
  std::vector<double> tensor_norms;

  double fission_source;
  double fission_source_prev_gen;

  const double err_k_tol;
  const double err_phi_tol;

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

  //TimerOutput                               computing_timer;

  std::vector<std_cxx11::shared_ptr<PreconditionSSOR<SparseMatrix<double> > > > pre_ho;
  SparsityPattern sparsity_pattern;
  ConstraintMatrix constraints;
};

#endif	// define  __ep_sn_h__

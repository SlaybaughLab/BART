#ifndef __problem_definition_h__
#define __problem_definition_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>

using namespace dealii;

/* Hard coded part, have to be changed in the future */
const unsigned int nmat = 50;
const unsigned int ngrp = 30;
const unsigned int z_levels = 30;
const unsigned int y_levels = 100;

template <int dim>
class ProblemDefinition
{
public:
  ProblemDefinition (ParameterHandler &prm);
  ~ProblemDefinition ();
  
  static void declare_parameters (ParameterHandler &prm);
  static std::string get_transport_model (ParameterHandler &prm);
  
  void process_input (ParameterHandler &prm);
  
  std::vector<std::vector<double> > get_sigma_t ();
  std::vector<std::vector<double> > get_inv_sigma_t ();
  std::vector<std::vector<double> > get_q ();
  std::vector<std::vector<double> > get_q_per_ster ();
  std::vector<std::vector<double> > get_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s_per_ster ();
  std::vector<std::vector<std::vector<double> > > get_ksi_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_ksi_nusigf_per_ster ();
  std::map<std::vector<unsigned int>, unsigned int> get_id_map ();
  std::unordered_map<unsigned int, bool> get_reflective_bc_map ();
  std::unordered_map<unsigned int, bool> get_fissile_id_map ();
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> get_component_index_map ();
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > get_inv_component_map ();
  
  std::vector<double> get_axis_maxes ();
  std::vector<unsigned int> get_ncells ();
  std::vector<double> get_cell_sizes ();
  
  unsigned int get_sn_order ();
  unsigned int get_n_dir ();
  unsigned int get_n_group ();
  unsigned int get_n_total_ho_vars ();
  std::vector<double> get_angular_weights ();
  std::vector<Tensor<1, dim> > get_all_directions ();
  std::vector<double> get_tensor_norms ();
  
  void initialize_component_index ();
  void print_angular_quad ();
  
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> component_index;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > inverse_component_index;
  
  bool get_nda_bool ();
  bool get_eigen_problem_bool ();
  bool get_reflective_bool ();
  bool get_print_sn_quad_bool ();
  unsigned int get_fe_order ();
  unsigned int get_uniform_refinement ();
  std::string get_discretization ();
  
  bool get_explicit_reflective_bool ();
  void initialize_ref_bc_index ();
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> get_reflective_direction_index_map ();
  
private:
  void process_problem_definition (ParameterHandler &prm);
  void preprocess_reflective_bc (ParameterHandler &prm);
  void process_coordinate_information (ParameterHandler &prm);
  void process_material_properties (ParameterHandler &prm);
  void process_eigen_material_properties (ParameterHandler &prm);
  void produce_angular_quad ();
  void initialize_relative_position_to_id_map (ParameterHandler &prm);
  
  double pi;
  
  unsigned int total_angle;
  
  std::string discretization;
  unsigned int n_azi;
  unsigned int n_group;
  unsigned int n_material;
  std::set<unsigned int> fissile_ids;
  
  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  unsigned int p_order;
  unsigned int global_refinements;
  
  std::unordered_map<unsigned int, bool> is_material_fissile;
  std::unordered_map<unsigned int, bool> is_reflective_bc;
  std::map<std::vector<unsigned int>, unsigned int> relative_position_to_id;
  
  bool do_print_sn_quad;
  bool is_explicit_reflective;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> reflective_direction_index;
  
  unsigned int n_dir;
  std::vector<double> wi;
  std::vector<Tensor<1, dim> > omega_i;
  std::vector<Tensor<1, 3> > omega_with_mu;
  std::vector<double> tensor_norms;
  unsigned int n_total_ho_vars;
  
  std::vector<unsigned int> ncell_per_dir;
  std::vector<double> cell_size_all_dir;
  std::vector<double> axis_max_values;
  
  std::vector<std::vector<double> > all_sigt;
  std::vector<std::vector<double> > all_inv_sigt;
  std::vector<std::vector<double> > all_ksi;
  std::vector<std::vector<double> > all_nusigf;
  std::vector<std::vector<double> > all_q;
  std::vector<std::vector<double> > all_q_per_ster;
  
  std::vector<std::vector<std::vector<double> > > all_sigs;
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf;
  std::vector<std::vector<std::vector<double> > > all_ksi_nusigf_per_ster;
  std::vector<std::vector<std::vector<double> > > ho_scaled_fiss_transfer_per_ster;
  std::vector<std::vector<std::vector<double> > > lo_scaled_fiss_transfer;
};


#endif	// define  __properties_h__

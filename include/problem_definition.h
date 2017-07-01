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

class ProblemDefinition
{
public:
  ProblemDefinition (ParameterHandler &prm);
  ~ProblemDefinition ();
  
  static void declare_parameters (ParameterHandler &prm);
  
  /*
   functions used to retrieve private members of ProblemDefinition<dim>
   */
  std::string get_transport_model ();// overloaded function
  std::string get_output_namebase ();
  std::string get_discretization ();
  std::string get_aq_name ();
  bool get_nda_bool ();
  bool get_eigen_problem_bool ();
  bool get_reflective_bool ();
  bool get_print_sn_quad_bool ();
  bool get_generated_mesh_bool ();
  bool get_explicit_reflective_bool ();
  unsigned int get_sn_order ();
  unsigned int get_n_dir ();
  unsigned int get_n_group ();
  unsigned int get_fe_order ();
  unsigned int get_uniform_refinement ();
  
private:
  std::string transport_model_name;
  std::string aq_name;
  std::string discretization;
  std::string mesh_filename;
  std::string output_namebase;
  bool is_mesh_generated;
  bool do_print_sn_quad;
  bool is_explicit_reflective;
  bool is_eigen_problem;
  bool do_nda;
  bool have_reflective_bc;
  unsigned int n_azi;
  unsigned int n_group;
  unsigned int n_material;
  unsigned int p_order;
  unsigned int global_refinements;
};


#endif	// define  __properties_h__

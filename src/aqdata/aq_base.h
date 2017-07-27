#ifndef __aq_base_h__
#define __aq_base_h__

#include <deal.II/base/tensor.h>
#include <deal.II/base/parameter_handler.h>

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

using namespace dealii;

template <int dim>
class AQBase
{
public:
  AQBase<dim> (ParameterHandler &prm);
  virtual ~AQBase ();

  virtual void produce_angular_quad ();
  virtual void initialize_component_index ();
  void print_angular_quad ();

  unsigned int get_sn_order ();
  unsigned int get_n_dir ();
  unsigned int get_n_total_ho_vars ();
  std::vector<double> get_angular_weights ();
  std::vector<double> get_tensor_norms ();
  std::vector<Tensor<1, dim> > get_all_directions ();
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  get_component_index_map ();
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
  get_inv_component_map ();
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  get_reflective_direction_index_map ();

protected:

  const double pi;
  double total_angle;
  bool have_reflective_bc;
  std::string transport_model_name;
  std::string discretization;
  unsigned int n_azi;
  unsigned int n_group;
  unsigned int n_dir;
  unsigned int n_total_ho_vars;
  std::vector<Tensor<1, dim> > omega_i;
  std::vector<double> wi;
  std::vector<double> tensor_norms;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  component_index;
  std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >
  inverse_component_index;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
  reflective_direction_index;

private:
  void make_aq (ParameterHandler &prm);
  std::string produce_quadrature_name ();
  void initialize_ref_bc_index ();

  std::string aq_name;
};

#endif //__aq_base_h__

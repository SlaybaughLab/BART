#ifndef __material_properties_h__
#define __material_properties_h__

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>

using namespace dealii;

class MaterialProperties
{
public:
  MaterialProperties (ParameterHandler &prm);
  ~MaterialProperties ();
  
  bool get_eigen_problem_bool ();
  unsigned int get_n_group ();
  unsigned int get_n_material ();
  
  std::unordered_map<unsigned int, bool> get_fissile_id_map ();
  
  std::vector<std::vector<double> > get_sigma_t ();
  std::vector<std::vector<double> > get_inv_sigma_t ();
  std::vector<std::vector<double> > get_q ();
  std::vector<std::vector<double> > get_q_per_ster ();
  std::vector<std::vector<double> > get_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s_per_ster ();
  std::vector<std::vector<std::vector<double> > > get_ksi_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_ksi_nusigf_per_ster ();
  
private:
  void process_material_properties (ParameterHandler &prm);
  void process_eigen_material_properties (ParameterHandler &prm);
  
  const double pi;
  
  bool is_eigen_problem;
  bool do_nda;
  unsigned int n_group;
  unsigned int n_material;
  
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::set<unsigned int> fissile_ids;
  
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
};

#endif //__material_properties_h__

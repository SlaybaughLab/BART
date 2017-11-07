#ifndef __material_properties_h__
#define __material_properties_h__

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

using namespace dealii;

//! This class read in and pre-process material properties.
/*!
 \author Weixiong Zheng
 \date 2017/04~06
 \todo Change scattering data from std::vector<std::vector<std::vector<double> > >
 to std::vector<Vector<double> > for scattering matrix eigenvalue decomposition.
 */
class MaterialProperties
{
public:
  /*!
   Class constructor.
   */
  MaterialProperties (ParameterHandler &prm);
  
  //! Class destructor.
  ~MaterialProperties ();
  
  unsigned int get_n_group ();
  unsigned int get_n_material ();
  
  /*!
   A function to retrieve mapping: material id->if material is fissile boolean.
   
   \return A Hash table representing the desirable mapping.
   */
  std::unordered_map<unsigned int, bool> get_fissile_id_map ();
  
  /*!
   A function to retrieve \f$\sigma_\mathrm{t}\f$ for all groups.
   */
  std::vector<std::vector<double> > get_sigma_t ();
  std::vector<std::vector<double> > get_inv_sigma_t ();
  std::vector<std::vector<double> > get_q ();
  std::vector<std::vector<double> > get_q_per_ster ();
  std::vector<std::vector<double> > get_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s ();
  std::vector<std::vector<std::vector<double> > > get_sigma_s_per_ster ();
  std::vector<std::vector<std::vector<double> > > get_chi_nusigf ();
  std::vector<std::vector<std::vector<double> > > get_chi_nusigf_per_ster ();
  
private:
  /*!
   Function to process all material properties. For eigenvalue problems, it will
   call this->process_eigen_material_properties to process fission-related
   properties.
   
   \param prm ParameterHandler object.
   \return Void.
   */
  void process_material_properties (ParameterHandler &prm);
  
  /*!
   Function to process eigenvalue problem specific material properties such as
   
   */
  void process_eigen_material_properties (ParameterHandler &prm);
  
  const double pi;//!< 3.14159...
  
  bool is_eigen_problem;//!< Boolean to determine if it's eigenvalue problem.
  bool do_nda;//!< Boolean to determine if NDA is used.
  unsigned int n_group;//!< Number of groups.
  unsigned int n_material;//!< Number of materials.
  
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::unordered_set<unsigned int> fissile_ids;//!< Set of fissile material IDs.
  
  std::vector<std::vector<double> > all_sigt;//!< \f$\sigma_\mathrm{t}\f$ for all groups.
  std::vector<std::vector<double> > all_inv_sigt;
  std::vector<std::vector<double> > all_chi;
  std::vector<std::vector<double> > all_nusigf;
  std::vector<std::vector<double> > all_q;
  std::vector<std::vector<double> > all_q_per_ster;
  
  std::vector<std::vector<std::vector<double> > > all_sigs;
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf;
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf_per_ster;
};

#endif //__material_properties_h__

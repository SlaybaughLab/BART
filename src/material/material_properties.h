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
 \todo Add functionality to perform eigenvalue decomposition.
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

  /*!
   A function to retrieve mapping: material id->if material is fissile boolean.
   
   \return A Hash table representing the desirable mapping.
   */
  std::unordered_map<unsigned int, bool> get_fissile_id_map ();
  
  //! A function to retrieve all \f$\sigma_\mathrm{t}\f$ for all groups.
  std::vector<std::vector<double> > get_sigma_t ();
  
  //! A function to retrieve all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  std::vector<std::vector<double> > get_inv_sigma_t ();
  
  //! A function to retrieve all fixed source value \f$Q\f$'s for all groups.
  std::vector<std::vector<double> > get_q ();
  
  //! A function to retrieve all \f$Q/(4\pi)\f$'s for all groups.
  std::vector<std::vector<double> > get_q_per_ster ();
  
  //! A function to retrieve all \f$\nu\sigma_\mathrm{f}\f$'s.
  std::vector<std::vector<double> > get_nusigf ();
  
  //! A function to retrieve all scattering transfer matrices.
  std::vector<std::vector<std::vector<double> > > get_sigma_s ();
  
  //! A function to retrieve all scattering transfer matrices scaled by \f$4\pi\f$.
  std::vector<std::vector<std::vector<double> > > get_sigma_s_per_ster ();
  
  //! A function to retrieve all \f$\chi\nu\sigma_\mathrm{f}\f$.
  std::vector<std::vector<std::vector<double> > > get_chi_nusigf ();
  
  //! A function to retrieve all \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$.
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
  
  //! Hash table to determine if a material is fissile.
  std::unordered_map<unsigned int, bool> is_material_fissile;
  
  std::unordered_set<unsigned int> fissile_ids;//!< Set of fissile material IDs.
  
  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::vector<std::vector<double> > all_sigt;
  
  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials
  std::vector<std::vector<double> > all_inv_sigt;
  
  /*!
   \f$\chi\f$ for all materials. A mistake when designing this is that it was 
   treated as other properties. Yet, it is group independent.
   
   \todo Change this to std::vector<double> chi, or using Hash table.
   */
  std::vector<std::vector<double> > all_chi;
  
  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all materials.
  std::vector<std::vector<double> > all_nusigf;
  
  //! Fixed source value \f$Q\f$'s of all groups for all materials.
  std::vector<std::vector<double> > all_q;
  
  //! Fixed source value \f$Q/(4\pi)\f$'s of all groups for all materials.
  std::vector<std::vector<double> > all_q_per_ster;
  
  //! Scattering transfer matrices for all materials.
  std::vector<std::vector<std::vector<double> > > all_sigs;
  
  //! Scattering transfer matrices scaled by \f$4\pi\f$
  std::vector<std::vector<std::vector<double> > > all_sigs_per_ster;
  
  //! \f$\chi\nu\sigma_\mathrm{f}\f$ for all materials.
  /*!
   It is pre-computed and transformed into a form of transfer matrix for saving
   computations.
   */
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf;
  
  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all materials.
  std::vector<std::vector<std::vector<double> > > all_chi_nusigf_per_ster;
};

#endif //__material_properties_h__

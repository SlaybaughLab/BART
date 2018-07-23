#ifndef BART_SRC_MATERIAL_MATERIALS_H_
#define BART_SRC_MATERIAL_MATERIALS_H_

#include <string>
#include <unordered_map>
#include <vector>
#include <unordered_set>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/full_matrix.h>

//! This class read in and pre-process material properties.
/*!
 \author Weixiong Zheng
 \date 2017/04~06

 \todo Change scattering data from std::vector<std::vector<std::vector<double> > >
 to std::vector<Vector<double> > for scattering matrix eigenvalue decomposition.
 \todo Add functionality to perform eigenvalue decomposition.
 */
class Materials {
 public:
  /*!
   Class constructor.
   */
  Materials (dealii::ParameterHandler &prm);

  //! Class destructor.
  ~Materials ();

  /*!
   A function to retrieve mapping: material id->if material is fissile boolean.

   \return A Hash table representing the desirable mapping.
   */
  std::unordered_map<int, bool> GetFissileIDMap () const;

  //! A function to retrieve all \f$\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetSigT () const;

  //! A function to retrieve all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  std::unordered_map<int, std::vector<double>> GetInvSigT () const;

  //! A function to retrieve all fixed source value \f$Q\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQ () const;

  //! A function to retrieve all \f$Q/(4\pi)\f$'s for all groups.
  std::unordered_map<int, std::vector<double>> GetQPerSter () const;

  //! A function to retrieve all \f$\nu\sigma_\mathrm{f}\f$'s.
  std::unordered_map<int, std::vector<double>> GetNuSigf () const;

  //! A function to retrieve all scattering transfer matrices.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigS () const;

  //! A function to retrieve all scattering transfer matrices scaled by \f$4\pi\f$.
  std::unordered_map<int, dealii::FullMatrix<double>> GetSigSPerSter () const;

  //! A function to retrieve all \f$\chi\nu\sigma_\mathrm{f}\f$.
  std::unordered_map<int, dealii::FullMatrix<double>> GetFissTransfer () const;

  //! A function to retrieve all \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$.
  std::unordered_map<int, dealii::FullMatrix<double>> GetFissTransferPerSter () const;

 private:
  /*!
   Function to process all material properties. For eigenvalue problems, it will
   call this->process_eigen_MATERIALS to process fission-related
   properties.

   \param prm ParameterHandler object.
   \return Void.
   */
  void ProcessMaterials (dealii::ParameterHandler &prm);

  /*!
   Function to process eigenvalue problem specific material properties such as

   */
  void ProcessEigenMaterials (dealii::ParameterHandler &prm);

  bool is_eigen_problem_;//!< Boolean to determine if it's eigenvalue problem.
  bool do_nda_;//!< Boolean to determine if NDA is used.
  int n_group_;//!< Number of groups.
  int n_material_;//!< Number of materials.

  //! Hash table to determine if a material is fissile.
  std::unordered_map<int, bool> is_material_fissile_;

  std::unordered_set<int> fissile_ids_;//!< Set of fissile material IDs.

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> sigt_;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials
  std::unordered_map<int, std::vector<double>> inv_sigt_;

  /*!
   \f$\chi\f$ for all materials. A mistake when designing this is that it was
   treated as other properties. Yet, it is group independent.

   \todo Change this to std::vector<double> chi, or using Hash table.
   */
  std::unordered_map<int, std::vector<double>> chi_;

  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> nu_sigf_;

  //! Fixed source value \f$Q\f$'s of all groups for all materials.
  std::unordered_map<int, std::vector<double>> q_;

  //! Fixed source value \f$Q/(4\pi)\f$'s of all groups for all materials.
  std::unordered_map<int, std::vector<double>> q_per_ster_;

  //! Scattering transfer matrices for all materials.
  std::unordered_map<int, dealii::FullMatrix<double>> sigs_;

  //! Scattering transfer matrices scaled by \f$4\pi\f$
  std::unordered_map<int, dealii::FullMatrix<double>> sigs_per_ster_;

  //! \f$\chi\nu\sigma_\mathrm{f}\f$ for all materials.
  /*!
   It is pre-computed and transformed into a form of transfer matrix for saving
   computations.
   */
  std::unordered_map<int, dealii::FullMatrix<double>> fiss_transfer_;

  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all materials.
  std::unordered_map<int, dealii::FullMatrix<double>> fiss_transfer_per_ster_;
};

#endif //BART_SRC_MATERIAL_MATERIALS_H_

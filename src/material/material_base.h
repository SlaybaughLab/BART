#ifndef BART_SRC_MATERIAL_MATERIAL_BASE_H_
#define BART_SRC_MATERIAL_MATERIAL_BASE_H_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

//! Defines the interface for Material classes providing material properties
/*!
  \author Joshua Rehak
  \date 2018/10/01
*/

class MaterialBase {
 public:

  virtual ~MaterialBase() = default;
  /*!
    returns an unordered_map from material ID to a
    boolean that is true if the material was labeled fissile
  */
  virtual std::unordered_map<int, bool> GetFissileIDMap() const = 0;

  //! Returns all \f$\sigma_\mathrm{t}\f$ for all groups.
  virtual std::unordered_map<int, std::vector<double>> GetSigT() const = 0;

  //! Returns all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  virtual std::unordered_map<int, std::vector<double>> GetInvSigT() const = 0;

  //! Returns all fixed source value \f$Q\f$'s for all groups.
  virtual std::unordered_map<int, std::vector<double>> GetQ() const = 0;

  //! Returns all \f$Q/(4\pi)\f$'s for all groups.
  virtual std::unordered_map<int, std::vector<double>> GetQPerSter() const = 0;

  //! Returns all \f$\nu\sigma_\mathrm{f}\f$'s.
  virtual std::unordered_map<int, std::vector<double>> GetNuSigF() const = 0;

  //! Returns all scattering transfer matrices.
  virtual std::unordered_map<int, dealii::FullMatrix<double>> GetSigS() const = 0;

  //! Returns all scattering transfer matrices scaled by \f$4\pi\f$.
  virtual std::unordered_map<int, dealii::FullMatrix<double>> GetSigSPerSter() const = 0;

  //! Returns \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  virtual std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigF() const = 0;

  //! Returns \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  virtual std::unordered_map<int, dealii::FullMatrix<double>> GetChiNuSigFPerSter() const = 0;

  /*!
   A function used to returns the eigen vector corresponding to the dominant mode
   of a dense matrix.

   \param dealii::FullMatrix the dense matrix that needs docomposing
   \return dealii::Vector A vector representing the eigenvector corresponding to
   the dominant mode (with largest-magnitude eigenvalue).

   \author Weixiong Zheng
   \date 2018/11

   \todo It is a naive implementation of the power iteration for arbitrary dense
   square matrix, given the fact that the scattering matrix will always have a
   dominant mode that is real. So the corresponding eigenvector has all entries
   to be real. A nice thing to have is to incorporate LAPACK in. deal.II has it
   to extract eigenvalues, but not eigenvectors. We either do a wrapping by ourselves
   with LAPACK, or put a ticket on deal.II github website and do a PR for it.
   */
  static dealii::Vector<double> GetEigenVectors (
      const dealii::FullMatrix<double>& mat) const;
};

#endif // BART_SRC_MATERIAL_MATERIAL_BASE_H_

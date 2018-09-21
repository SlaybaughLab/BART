#ifndef BART_SRC_MATERIAL_MATERIAL_PROPERTIES_I_H_
#define BART_SRC_MATERIAL_MATERIAL_PROPERTIES_I_H_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

class MaterialPropertiesI {
 public:

  virtual ~MaterialPropertiesI() = default;
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
  
};

#endif

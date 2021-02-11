#ifndef BART_SRC_MATERIAL_MATERIAL_I_HPP_
#define BART_SRC_MATERIAL_MATERIAL_I_HPP_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

//! Material properties
namespace bart::material {

//! Defines the interface for Material classes providing material properties
/*!
  \author Joshua Rehak
  \date 2018/10/01
*/
class MaterialI {
 public:
  //! Dealii full matrix type
  using DealiiMatrix = dealii::FullMatrix<double>;
  //! Maps material ID (int) to another type
  template <typename MappedType>
  using MaterialIDMappedTo = std::unordered_map<int, MappedType>;
  virtual ~MaterialI() = default;
  /*!
    returns an unordered_map from material ID to a
    boolean that is true if the material was labeled fissile
  */
  virtual auto GetFissileIDMap() const -> MaterialIDMappedTo<bool> = 0;
  //! Returns diffusion coefficient or 0 if not provided
  virtual auto GetDiffusionCoef() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Returns all \f$\sigma_\mathrm{t}\f$ for all groups.
  virtual auto GetSigT() const -> MaterialIDMappedTo<std::vector<double>>  = 0;
  //! Returns all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  virtual auto GetInvSigT() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Returns all fixed source value \f$Q\f$'s for all groups.
  virtual auto GetQ() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Returns all \f$Q/(4\pi)\f$'s for all groups.
  virtual auto GetQPerSter() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Returns all \f$\nu\sigma_\mathrm{f}\f$'s.
  virtual auto GetNuSigF() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Returns all scattering transfer matrices.
  virtual auto GetSigS() const -> MaterialIDMappedTo<DealiiMatrix> = 0;
  //! Returns all scattering transfer matrices scaled by \f$4\pi\f$.
  virtual auto GetSigSPerSter() const -> MaterialIDMappedTo<DealiiMatrix> = 0;
  //! Returns \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  virtual auto GetChiNuSigF() const -> MaterialIDMappedTo<DealiiMatrix> = 0;
  //! Returns \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  virtual auto GetChiNuSigFPerSter() const -> MaterialIDMappedTo<DealiiMatrix> = 0;
};

} // namespace bart::material

#endif // BART_SRC_MATERIAL_MATERIAL_I_HPP_

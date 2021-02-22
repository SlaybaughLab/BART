#ifndef BART_SRC_DATA_CROSS_SECTIONS_CROSS_SECTIONS_I_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_CROSS_SECTIONS_I_HPP_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

//! \brief Nuclear cross-section data
namespace bart::data::cross_sections {

class CrossSectionsI {
 public:
  using DealiiMatrix = dealii::FullMatrix<double>;
  //! Maps material ID (int) to another type
  template <typename MappedType>
  using MaterialIDMappedTo = std::unordered_map<int, MappedType>;

  virtual ~CrossSectionsI() = default;

  //! Diffusion coefficient of all groups for all materials.
  virtual auto diffusion_coef() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  virtual auto sigma_t() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  virtual auto inverse_sigma_t() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! Scattering matrices for all materials (i.e. \f$\sigma_\mathrm{s,g'\to g}\f$).
  virtual auto sigma_s() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> = 0;
  //! \f$\sigma_\mathrm{s,g'\to g}/(4\pi)\f$
  virtual auto sigma_s_per_ster() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> = 0;
  //! \f$Q\f$ values of all groups for all materials.
  virtual auto q() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! \f$Q/(4\pi)\f$ values of all groups for all materials.
  virtual auto q_per_ster() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  virtual auto is_material_fissile() const -> MaterialIDMappedTo<bool> = 0;
  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all fissile materials.
  virtual auto nu_sigma_f() const -> MaterialIDMappedTo<std::vector<double>> = 0;
  //! \f$\chi\nu\sigma_\mathrm{f}\f$ of all incident and outgoing groups for fissile materials.
  virtual auto fiss_transfer() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> = 0;
  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for fissile materials.
  virtual auto fiss_transfer_per_ster() const -> MaterialIDMappedTo<dealii::FullMatrix<double>> = 0;
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_CROSS_SECTIONS_I_HPP_

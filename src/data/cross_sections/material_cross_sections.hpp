#ifndef BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_

#include "cross_sections_i.hpp"

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "data/material/material_i.hpp"

namespace bart::data::cross_sections {

/*! \brief Cross-sections generated from a material.
 *
 * This is the default implementation of cross-sections that parses a MaterialI object for data.
 *
 */
class MaterialCrossSections: public CrossSectionsI {
 public:
  using CrossSectionsI::DealiiMatrix;
  using CrossSectionsI::MaterialIDMappedTo;
  
  MaterialCrossSections(data::material::MaterialI &materials);

  auto diffusion_coef() const -> MaterialIDMappedTo<std::vector<double>> override { return diffusion_coef_; }
  auto sigma_t() const -> MaterialIDMappedTo<std::vector<double>> override { return sigma_t_; }
  auto inverse_sigma_t() const -> MaterialIDMappedTo<std::vector<double>> override { return inverse_sigma_t_; }
  auto sigma_s() const -> MaterialIDMappedTo<DealiiMatrix> override { return sigma_s_; }
  auto sigma_s_per_ster() const -> MaterialIDMappedTo<DealiiMatrix> override { return sigma_s_per_ster_; }
  auto q() const -> MaterialIDMappedTo<std::vector<double>> override { return q_; }
  auto q_per_ster() const -> MaterialIDMappedTo<std::vector<double>> override { return q_per_ster_; }
  auto is_material_fissile() const -> MaterialIDMappedTo<bool> override { return is_material_fissile_; }
  auto nu_sigma_f() const -> MaterialIDMappedTo<std::vector<double>> override { return nu_sigma_f_; }
  auto fiss_transfer() const -> MaterialIDMappedTo<DealiiMatrix> override { return fiss_transfer_; }
  auto fiss_transfer_per_ster() const -> MaterialIDMappedTo<DealiiMatrix> override { return fiss_transfer_per_ster_; }
 private:
  const MaterialIDMappedTo<std::vector<double>> diffusion_coef_;
  const MaterialIDMappedTo<std::vector<double>> sigma_t_;
  const MaterialIDMappedTo<std::vector<double>> inverse_sigma_t_;
  const MaterialIDMappedTo<DealiiMatrix> sigma_s_;
  const MaterialIDMappedTo<DealiiMatrix> sigma_s_per_ster_;
  const MaterialIDMappedTo<std::vector<double>> q_;
  const MaterialIDMappedTo<std::vector<double>> q_per_ster_;
  const MaterialIDMappedTo<bool> is_material_fissile_;
  const MaterialIDMappedTo<std::vector<double>> nu_sigma_f_;
  const MaterialIDMappedTo<DealiiMatrix> fiss_transfer_;
  const MaterialIDMappedTo<DealiiMatrix> fiss_transfer_per_ster_;
};

} // namespace bart::data::cross_sections

#endif // BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_

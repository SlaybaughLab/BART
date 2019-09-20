#ifndef BART_SRC_DATA_CROSS_SECTIONS_H_
#define BART_SRC_DATA_CROSS_SECTIONS_H_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "../material/material_base.h"

namespace bart {

namespace data {

struct CrossSections {
  typedef int MaterialID;
  
  CrossSections(MaterialBase &materials);

  //! Diffusion coefficient of all groups for all materials.
  const std::unordered_map<MaterialID, std::vector<double>> diffusion_coef;

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<MaterialID, std::vector<double>> sigma_t;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<MaterialID, std::vector<double>> inverse_sigma_t;

  //! Scattering matrices for all materials (i.e. \f$\sigma_\mathrm{s,g'\to g}\f$).
  const std::unordered_map<MaterialID, dealii::FullMatrix<double>> sigma_s;

  //! \f$\sigma_\mathrm{s,g'\to g}/(4\pi)\f$
  const std::unordered_map<MaterialID, dealii::FullMatrix<double>> sigma_s_per_ster;

  //! \f$Q\f$ values of all groups for all materials.
  const std::unordered_map<MaterialID, std::vector<double>> q;

  //! \f$Q/(4\pi)\f$ values of all groups for all materials.
  const std::unordered_map<MaterialID, std::vector<double>> q_per_ster;

  const std::unordered_map<MaterialID, bool> is_material_fissile;

  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all fissile materials.
  const std::unordered_map<MaterialID, std::vector<double>> nu_sigma_f;

  //! \f$\chi\nu\sigma_\mathrm{f}\f$ of all incident and outgoing groups for fissile materials.
  const std::unordered_map<MaterialID, dealii::FullMatrix<double>> fiss_transfer;

  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for fissile materials.
  const std::unordered_map<MaterialID, dealii::FullMatrix<double>> fiss_transfer_per_ster;
  
}; 
  
} // namespace data

} // namespace bart 

#endif // BART_SRC_DATA_CROSS_SECTIONS_H_

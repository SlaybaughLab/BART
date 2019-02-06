#ifndef BART_SRC_DATA_CROSS_SECTIONS_H_
#define BART_SRC_DATA_CROSS_SECTIONS_H_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

namespace bart {

namespace data {

struct CrossSections {
  virtual ~CrossSections() = default;

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> sigma_t;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> inverse_sigma_t;

  //! Scattering matrices for all materials (i.e. \f$\sigma_\mathrm{s,g'\to g}\f$).
  const std::unordered_map<int, dealii::FullMatrix<double>> sigma_s;

  //! \f$\sigma_\mathrm{s,g'\to g}/(4\pi)\f$
  const std::unordered_map<int, dealii::FullMatrix<double>> sigma_s_per_ster;
}; 
  
} // namespace data

} // namespace bart

#endif // BART_SRC_DATA_CROSS_SECTIONS_I_H_

#ifndef BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_
#define BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_

#include <memory>

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "calculator/cell/integrated_fission_source_i.hpp"
#include "domain/domain_types.h"
#include "utility/has_dependencies.h"
#include "utility/has_description.h"

namespace bart {

namespace domain {
template <int dim> class DefinitionI;
} // namespace domain

namespace calculator::cell {

/*! \brief Default implementation of the aggregated fission source using Dealii.
 *
 * This implementation of TotalAggregatedFissionSourceI uses the deal.II domain and an IntegratedFissionSource class
 * to calculate the aggregated fission source. The returned value is given by,
 *
 * \f[
 *   F_{\mathcal{D}} = \sum_{k = 0}^{K} F_{k}
 * \f]
 *
 * where \f$K\f$ is the total number of cells in the domain \f$\mathcal{D}\f$ and \f$F_{k}\f$ is the integrated
 * fission source for that cell.
 *
 * @tparam dim spatial dimension
 */
template <int dim>
class TotalAggregatedFissionSource : public TotalAggregatedFissionSourceI, public utility::HasDescription,
                                     public utility::HasDependencies {
 public:
  using Domain = domain::DefinitionI<dim>;
  using IntegratedFissionSource = IntegratedFissionSourceI<dim>;
  using SystemMoments = system::moments::SphericalHarmonicI;

  /*! \brief Constructor.
   *
   * @param cell_fission_source_ptr integrated fission source calculator for cells.
   * @param domain_ptr domain that contains the cells of interest for aggregating the fission source.
   */
  TotalAggregatedFissionSource(std::unique_ptr<IntegratedFissionSource> cell_fission_source_ptr,
                               std::shared_ptr<Domain> domain_ptr);
  virtual ~TotalAggregatedFissionSource() = default;

  [[nodiscard]] auto AggregatedFissionSource(SystemMoments*) const -> double override;

  /*! \brief Access integrated fission source calculator dependency. */
  IntegratedFissionSource* cell_fission_source_ptr() const { return cell_fission_source_ptr_.get(); }
   /*! \brief Access domain dependency. */
  Domain* domain_ptr() const { return domain_ptr_.get(); }
 private:
  std::unique_ptr<IntegratedFissionSource> cell_fission_source_ptr_;
  std::shared_ptr<Domain> domain_ptr_;
};

} // namespace calculator::cell

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TOTAL_AGGREGATED_FISSION_SOURCE_HPP_
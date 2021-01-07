#ifndef BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_I_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_I_HPP_

#include <map>

#include <deal.II/lac/vector.h>

#include "quadrature/quadrature_types.h"
#include "utility/named_type.h"

namespace bart::quadrature::calculators {
/*! \brief Interface for a class that integrates angular flux using a quadrature set.
 *
 * The provided functions calculate the net current at a degree of freedom \f$i\f$:
 * \f[
 * \vec{J}_i = \int \hat{\Omega} \psi(\hat{\Omega})_i d\hat{\Omega} = \sum_{m = 0}^M w_m\hat{\Omega}_m\psi_{i,m}\;,
 * \f]
 * the magnitude of the current in direction \f$\hat{n}\f$:
 * \f[
 * j_{\hat{n}} = \int_{\hat{n} \cdot \hat{\Omega} \ge 0} |\hat{n} \cdot \hat{\Omega}| \psi(\hat{\Omega})_i d\hat{\Omega} = \sum_{\Omega_m \mid \hat{n} \cdot \Omega \ge 0} w_m|\hat{n} \cdot \hat{\Omega}_m|\psi_{i,m}\;,
 * \f]
 * and the integrated angular flux in direction \f$\hat{n}\f$:
 *
 * \f[
 * \phi_{\hat{n}} = \int_{\hat{n} \cdot \hat{\Omega} \ge 0} \psi(\hat{\Omega})_i d\hat{\Omega} = \sum_{\Omega_m \mid \hat{n} \cdot \Omega \ge 0} w_m\psi_{i,m}\;.
 * \f]
 */
class AngularFluxIntegratorI {
 public:
  virtual ~AngularFluxIntegratorI() = default;
  using Vector = dealii::Vector<double>;
  using VectorPtr = std::shared_ptr<dealii::Vector<double>>;
  using VectorMap = std::map<quadrature::QuadraturePointIndex, VectorPtr>;
  using DegreeOfFreedom = utility::NamedType<int, struct DegreeOfFreedomParam>;

  virtual auto NetCurrent(const VectorMap&, const DegreeOfFreedom) const -> Vector = 0;
  virtual auto DirectionalCurrent(const VectorMap&, const Vector normal, const DegreeOfFreedom) const -> double = 0;
  virtual auto DirectionalFlux(const VectorMap&, const Vector normal, const DegreeOfFreedom) const -> double = 0;
};

} // namespace bart::quadrature::calculators

#endif //BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_I_HPP_
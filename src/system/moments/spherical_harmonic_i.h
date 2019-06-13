#ifndef BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_
#define BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_

#include "system/moments/spherical_harmonic_types.h"

namespace bart {

namespace system {

namespace moments {

/*! \brief Interface for spherical harmonic moments storage class
 *
 * This class is the interface for classes that stores system spherical harmonic
 * moments. Each is stored using an index made up of three values that corrspond to
 * group, degree, and order of the moment. That is, each is given by:
 * \f[
 *
 * \phi_{g}^{\ell, m}(\vec{r}) = \sum_{i = 0}^{N}w_i\psi_{g}(\vec{r}, \hat{\Omega}_i)Y_{\ell}^{m}(\hat{\Omega}_i)
 *
 * \f]
 *
 * where \f$N\f$ is the total number of angles in the angular quadrature, and
 * \f$Y_{\ell}^{m}\f$ are the spherical harmonics of degree
 * \f$ \{\ell \in \mathbb{Z} \mid 0 \leq \ell \leq \ell_{\text{max}}\} \f$
 * and order
 * \f$ \{m \in \mathbb{Z} \mid  |m| \leq \ell}\} \f$.
 */
class SphericalHarmonicI {
 public:
  virtual ~SphericalHarmonicI() = default;

  /*! \brief Returns the full mapping of moments.
   *
   * By default this interface does not provide a method for modifying the
   * map (it is returned as a constant).
   *
   */
  virtual const MomentsMap& moments() const = 0;

  virtual const MomentVector& operator[](const MomentIndex) const = 0;
  virtual MomentVector& operator[](const MomentIndex) = 0;


  /*! \brief Returns the total number of energy groups. */
  virtual int total_groups() const = 0;
  /*! \brief Returns the value of \f$\ell_{\text{max}}\f$. */
  virtual int max_harmonic_l() const = 0;
};

} // namespace moments

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_MOMENTS_SPHERICAL_HARMONIC_I_H_
#ifndef BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_
#define BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_

#include "quadrature/quadrature_types.h"

#include <utility>
#include <vector>

#include "utility/has_description.h"

namespace bart {

namespace quadrature {

/*! \brief Generates a quadrature set, generally for one quadrant.
 *
 * Quadrature points are generated using the quadrature types that make them
 * ideal for instantiating Ordinates and Quadrature Points.
 *
 * @tparam dim spatial dimension.
 */
template <int dim>
class QuadratureGeneratorI : public bart::utility::HasDescription {
 public:
  virtual ~QuadratureGeneratorI() = default;
  /*! \brief Generates a quadrature set.
   * Sets are generated as a vector containing pairs of cartesian positions and
   * weights. This is appropriate for the input to the
   * quadrature::utility::FillQuadratureSet function.
   */
  virtual std::vector<std::pair<CartesianPosition<dim>, Weight>>
  GenerateSet() const = 0;
  /*! \brief Returns the order of the quadrature generated.
   * The meaning of the order may vary by the specific type of quadrature.
   * @return int indicating the order.
   */
  virtual int order() const = 0;
};

} // namespace quadrature

} // namespace bart




#endif //BART_SRC_QUADRATURE_ANGULAR_ANGULAR_QUADRATURE_GENERATOR_I_H_

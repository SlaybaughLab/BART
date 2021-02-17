#ifndef BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_I_HPP_

#include <deal.II/lac/vector.h>

//! \brief Calculators for residuals
namespace bart::calculator::residual {

/*! \brief Interface for a residual calculator that uses the difference between two dealii vectors */
class VectorDifferenceI {
 public:
  using Vector = dealii::Vector<double>;
  virtual ~VectorDifferenceI() = default;
  /*! \brief Calculate the residual using two vectors.
   *  The residual is simply given by \f$\vec{R} = \vec{v}_m - \vec{v}_s\f$ where the subscripts \f$m\f$ and \f$s\f$
   *  indicate the minuend and subtrahend. The minuend is the first vector provided to this function, the
   *  subtrahend is the second.
   * @return residual vector.
   */
  virtual auto CalculateResidual(const Vector&, const Vector&) const -> Vector = 0;
  /*! \brief Calculate the residual and then weigh by some value.
   * The residual is calculated in a similar manner when called without weight
   * \f$\vec{R} = w\left(\vec{v}_m - \vec{v}_s\right)\f$ where \f$w\f$ is the weight.
   *
   * @param weight parameter to weigh the residual
   * @return weighed residual vector
   */
  virtual auto CalculateResidual(const Vector&, const Vector&, const double weight) const -> Vector = 0;
};


} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_VECTOR_DIFFERENCE_I_HPP_

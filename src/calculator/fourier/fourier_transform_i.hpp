#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

#include "utility/named_type.h"
#include "utility/has_description.h"

//! \brief Calculators for calculating the Fourier transform.
namespace bart::calculator::fourier {

using Normalized = utility::NamedType<bool, struct NormalizedStruct>;

/*! \brief Interface for classes that calculate the discrete fourier transform. */
class FourierTransformI : public utility::HasDescription {
 public:
  /// Vector of complex values, used for both input and output of the DFT.
  using ComplexVector = std::vector<std::complex<double>>;
  virtual ~FourierTransformI() = default;
  /*! \brief Calculate the DFT of a provided complex vector.
   *
   * @param input input vector of complex values to take the DFT of.
   * @param normalize_return if true, the result is normalized by dividing by the total size of the vector.
   * @return vector of complex vector that holds the DFT of the input.
   */
  virtual auto CalculateDFT(const ComplexVector& input, Normalized normalize_return) -> ComplexVector = 0;
  /*! \brief Calculate the DFT of a provided dealii-format Vector.
   *
   * @param input input vector of real values to take the DFT of.
   * @param normalize_return if true, the result is normalized by dividing by the total size of the vector.
   * @return vector of complex vector that holds the DFT of the input.
   */
  virtual auto CalculateDFT(const dealii::Vector<double>& input, Normalized normalize_return) -> ComplexVector = 0;
};

} // namespace bart::calculator::fourier

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_

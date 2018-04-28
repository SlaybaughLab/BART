#ifndef BART_SRC_COMMON_NUMBERS_H_
#define BART_SRC_COMMON_NUMBERS_H_

#include <deal.II/base/numbers.h>

//! This namespace provide constant numbers
/*!
 \author Weixiong Zheng
 \date 2018/04
 */
namespace bconst {
  static const kPi = dealii::numbers::PI;//!< Pi
  static const kTwoPi = 2.*dealii::numbers::PI;//!< 2Pi
  static const kFourPi = 4.*dealii::numbers::PI;//!< 4Pi
  static const kInvPi = 1./dealii::numbers::PI;//!< 1/Pi
  static const kInvTwoPi = 1./kTwoPi;//!< 1/(2*Pi)
  static const kInvFourPi = 1./kFourPi;//!< 1/(4*Pi)
  static const kSmall = 1.0e-15;//!< Non-machine-error small number.
}

#endif //BART_SRC_COMMON_NUMBERS_H_

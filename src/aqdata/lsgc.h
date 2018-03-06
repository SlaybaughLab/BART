#ifndef BART_SRC_AQDATA_LSGC_H__
#define BART_SRC_AQDATA_LSGC_H__

#include "aq_base.h"

//! This class produces level-symmetric Gauss-Chebyshev quadrature.
/*!
 This class is derived from AQBase<dim> with LSGC rule. The polar angular points
 are determined by 1D Gauss-Legendre quadrature. Azimuthally, the points are
 placed uniformly per polar level. Detailed descriptions can be found in literature
 such as <a href="http://oaktrust.library.tamu.edu/bitstream/handle/1969.1/ETD-T
 AMU-2010-12-8586/JARRELL-DISSERTATION.pdf?sequence=2" style="color:blue"><b>
 J. J. Jarrell's dissertation</b></a>.

 \author Weixiong Zheng
 \date 2017/04
 */
template <int dim>
class LSGC : public AQBase<dim> {
 public:
  /*!
   Class constructor.

   \param prm ParameterHandler object.
   */
  LSGC (dealii::ParameterHandler &prm);

  //!< Destructor.
  ~LSGC ();

  /*!
   This function overrides AQBase<dim>::produce_angular_quad () to specifically
   produce level-symmetric Gauss-Chebyshev quadrature.

   \return Void.
   */
  void ProduceAQ ();
};

#endif// BART_SRC_AQDATA_LSGC_H__

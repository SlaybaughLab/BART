#ifndef __lsgc_h__
#define __lsgc_h__

#include "aq_base.h"

using namespace dealii;

//! This class produces level-symmetric Gauss-Chebyshev quadrature.
template <int dim>
class LSGC : public AQBase<dim>
{
public:
  /*!
   Class constructor.
   
   \param prm ParameterHandler object.
   */
  LSGC (ParameterHandler &prm);
  
  //!< Destructor.
  ~LSGC ();
  
  /*!
   This function override AQBase<dim>::produce_angular_quad () to specifically
   produce level-symmetric Gauss-Chebyshev quadrature.
   
   \return Void.
   */
  void produce_angular_quad ();
};

#endif//__lsgc_h__

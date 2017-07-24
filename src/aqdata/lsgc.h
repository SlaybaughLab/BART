#ifndef __lsgc_h__
#define __lsgc_h__

#include "aq_base.h"

using namespace dealii;

template <int dim>
class LSGC : public AQBase<dim>
{
public:
  LSGC (ParameterHandler &prm);
  ~LSGC ();
  
  void produce_angular_quad ();
};

#endif//__aq_lsgc_h__

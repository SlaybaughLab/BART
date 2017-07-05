#ifndef __aq_lsgc_h__
#define __aq_lsgc_h__

#include "../base/aq_base.h"

using namespace dealii;

template <int dim>
class AQLSGC : public AQBase<dim>
{
public:
  AQLSGC (ParameterHandler &prm);
  ~AQLSGC ();
  
  void produce_angular_quad ();
};

#endif//__aq_lsgc_h__

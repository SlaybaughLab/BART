#ifndef __MODEL_MANAGER__H__
#define __MODEL_MANAGER__H__

#include "transport_base.h"
#include "even_parity.h"

template <int dim>
class ModelManager
{
public:
  ModelManager ();
  ~ModelManager ();

  static std_cxx11::shared_ptr<TransportBase<dim> >
  build_transport_model (std::string &transport_model_name,
                         ParameterHandler &prm);

};

#endif //__MODEL_MANAGER__H__

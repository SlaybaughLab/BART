#include "../../include/common/model_manager.h"
#include "../../include/transport/base/transport_base.h"
#include "../../include/transport/derived/even_parity.h"


ModelManager::ModelManager (ParameterHandler &prm)
:
transport_model_name(prm.get("transport model")),
dim(prm.get_integer("problem dimension"))
{
  // register new methods here
  method_index["ep"] = 0;
}

ModelManager::~ModelManager ()
{
}

void ModelManager::build_and_run_model (ParameterHandler &prm)
{
  AssertThrow(dim!=1,
              ExcMessage("1D is not implemented"));
  
  unsigned int mi = method_index[transport_model_name];
  switch (mi)
  {
    case 0:
    {
      if (dim==2)
      {
        std_cxx11::shared_ptr<TransportBase<2> > tb = std_cxx11::shared_ptr<TransportBase<2> > (new EvenParity<2>(prm));
        tb->run ();
      }
      else
      {
        std_cxx11::shared_ptr<TransportBase<3> > tb = std_cxx11::shared_ptr<TransportBase<3> > (new EvenParity<3>(prm));
        tb->run ();
      }
    }
      break;
      
    default:
      break;
  }
}

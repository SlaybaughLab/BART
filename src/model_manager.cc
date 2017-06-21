#include "../include/model_manager.h"

template <int dim>
ModelManager<dim>::ModelManager ()
{
}

template <int dim>
ModelManager<dim>::~ModelManager ()
{
}

template <int dim>
std_cxx11::shared_ptr<TransportBase<dim> >
ModelManager<dim>::build_transport_model (std::string &transport_model_name,
                       ParameterHandler &prm)
{
  if (transport_model_name=="ep")
    return std_cxx11::shared_ptr<TransportBase<dim> > (new EvenParity<dim>(prm));
}

template class ModelManager<2>;
template class ModelManager<3>;

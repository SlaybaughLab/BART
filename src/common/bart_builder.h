#ifndef __bart_builder_h__
#define __bart_builder_h__

#include <deal.II/base/parameter_handler.h>

#include <string>
#include <map>

#include "../transport/transport_base.h"

using namespace dealii;

template <int dim>
std_cxx11::shared_ptr<TransportBase<dim> >
build_transport_model (ParameterHandler &prm);

template <int dim>
std_cxx11::shared_ptr<AQBase<dim> >
build_aq_model (ParameterHandler &prm);

#endif //__bart_builder_h__

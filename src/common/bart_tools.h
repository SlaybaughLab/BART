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

void radio (std::string str);

void radio (std::string str1, std::string str2);

void radio (std::string str, double num);

void radio (std::string str1, unsigned int num1,
            std::string str2, unsigned int num2,
            std::string str3, unsigned int num3);

void radio (std::string str, unsigned int num);

void radio (std::string str, bool boolean);

void radio ();

#endif //__bart_builder_h__

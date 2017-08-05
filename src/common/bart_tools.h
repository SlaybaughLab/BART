#ifndef __bart_builder_h__
#define __bart_builder_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe_poly.h>

#include <string>
#include <map>

#include "../equations/equation_base.h"

using namespace dealii;

template <int dim>
FE_Poly<TensorProductPolynomials<dim>,dim,dim>*
build_finite_element (ParameterHandler &prm);

template <int dim>
std_cxx11::shared_ptr<MeshGenerator<dim> >
build_mesh (ParameterHandler &prm);

std_cxx11::shared_ptr<MaterialProperties> build_material (ParameterHandler &prm);

template <int dim>
std_cxx11::shared_ptr<IterationBase<dim> >
build_iterative_solver (ParameterHandler &prm,
                        const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                        const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr);

template <int dim>
std_cxx11::shared_ptr<TransportBase<dim> >
build_transport_model (ParameterHandler &prm,
                       const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                       const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                       const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);

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

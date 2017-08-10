#ifndef __bart_tools_h__
#define __bart_tools_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/fe/fe_poly.h>

#include <string>
#include <map>

#include "../equations/equation_base.h"
#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"
#include "../common/material_properties.h"

using namespace dealii;

/** \brief Function to build finite element for general dimensions specified by
 * user.
 *
 * \parameter prm A reference to processed ParameterHandler instance
 * \return A raw pointer of finite elements derived from FE_Poly
 */
template <int dim>
FE_Poly<TensorProductPolynomials<dim>,dim,dim>*
build_finite_element (ParameterHandler &prm);

/** \brief Function to build mesh in calculations for general dimensions
 *
 * \parameter prm A reference to processed ParameterHandler instance
 * \return A shared pointer to MeshGenerator<dim> instance
 */
template <int dim>
std_cxx11::shared_ptr<MeshGenerator<dim> >
build_mesh (ParameterHandler &prm);

/** \brief Function to build pointer to MaterialProperties class.
 *
 * \parameter prm A reference to processed ParameterHandler instance
 * \return A shared pointer to MaterialProperties instance
 */
std_cxx11::shared_ptr<MaterialProperties> build_material (ParameterHandler &prm);

/** \brief Build specific transport model
 *
 * \parameter msh_ptr shared_ptr for MeshGenerator<dim> instance
 * \parameter aqd_ptr shared_ptr for AQBase<dim> instance
 * \parameter mat_ptr shared_ptr for MaterialProperties instance
 * \return a shared pointer to transport model
 */
template <int dim>
std_cxx11::shared_ptr<EquationBase<dim> >
build_transport_model (ParameterHandler &prm,
                       const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
                       const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
                       const std_cxx11::shared_ptr<MaterialProperties> mat_ptr);

/** \brief Function to build angular quadrature for general dimensions
 *
 * \parameter prm A processed ParameterHandler instance
 * \return a shared pointer to specific angular quadrature
 */
template <int dim>
std_cxx11::shared_ptr<AQBase<dim> >
build_aq_model (ParameterHandler &prm);

/** \brief ConditionalOStream object used to output things on screen with processor 0
 */
ConditionalOStream pout (std::cout,
                         Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0);

void radio (std::string str);

void radio (std::string str1, std::string str2);

void radio (std::string str, double num);

void radio (std::string str1, unsigned int num1,
            std::string str2, unsigned int num2,
            std::string str3, unsigned int num3);

void radio (std::string str, unsigned int num);

void radio (std::string str, bool boolean);

void radio ();

#endif //__bart_tools_h__

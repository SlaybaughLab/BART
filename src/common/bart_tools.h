#ifndef __bart_tools_h__
#define __bart_tools_h__

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

#include <string>

#include "../equation/equation_base.h"
#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"
#include "../iteration/iterations.h"
#include "../iteration/eigen_base.h"
#include "../iteration/mg_base.h"
#include "../iteration/ig_base.h"
#include "../material/material_properties.h"
#include "../iteration/power_iteration.h"
#include "../iteration/gauss_sidel.h"
#include "../equation/even_parity.h"
#include "../aqdata/lsgc.h"

using namespace dealii;

template <int dim>
void build_finite_element
(FE_Poly<TensorProductPolynomials<dim>,dim,dim>* fe, ParameterHandler &prm)
{
  std::string discretization = prm.get ("spatial discretization");
  AssertThrow (discretization!="none",
               ExcMessage("valid spatial discretizations are dfem and cfem"))
  int p_order = prm.get_integer("finite element polynomial degree");
  if (discretization=="dfem")
    fe = new FE_DGQ<dim> (p_order);
  else
    fe = new FE_Q<dim> (p_order);
}

void build_linalg
(std_cxx11::shared_ptr<PreconditionerSolver> alg_ptr,
 const ParameterHandler &prm,
 std::string equation_name,
 unsigned int& n_total_vars)
{
  alg_ptr = std_cxx11::shared_ptr<PreconditionerSolver>
  (new PreconditionerSolver(prm, equation_name, n_total_vars));
}

template <int dim>
void build_iterations
(std_cxx11::shared_ptr<Iterations<dim> > itr_ptr,
 const ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
{
  itr_ptr = std_cxx11::shared_ptr<Iterations<dim> >
  (new Iterations<dim> (prm, msh_ptr, aqd_ptr, mat_ptr));
}

template <int dim>
void build_mesh
(std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr, ParameterHandler &prm)
{
  msh_ptr =
  std_cxx11::shared_ptr<MeshGenerator<dim> >(new MeshGenerator<dim>(prm));
}

void build_material
(std_cxx11::shared_ptr<MaterialProperties> mat_ptr, ParameterHandler &prm)
{
  mat_ptr =
  std_cxx11::shared_ptr<MaterialProperties> (new MaterialProperties (prm));
}

template <int dim>
void build_equation
(std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr,
 std::string equation_name,
 const ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
{
  // TODO: add NDA to it after having NDA class
  std_cxx11::shared_ptr<EquationBase<dim> > equation_pointer;
  if (equation_name=="ep")
    equ_ptr =
    std_cxx11::shared_ptr<EquationBase<dim> >
    (new EvenParity<dim> (equation_name, prm, msh_ptr, aqd_ptr, mat_ptr));
}

template <int dim>
void build_aq_model
(std_cxx11::shared_ptr<AQBase<dim> > aq_pointer, ParameterHandler &prm)
{
  std::string aq_name = prm.get ("angular quadrature name");
  AssertThrow (aq_name!="none",
               ExcMessage("angular quadrature name incorrect or missing"));
  //std_cxx11::shared_ptr<AQBase<dim> > aq_pointer;
  if (aq_name=="lsgc")
    aq_pointer = std_cxx11::shared_ptr<AQBase<dim> > (new LSGC<dim>(prm));
  //return aq_pointer;
}

/** \brief Function used to build pointer to instance of InGroupBase's derived class
 */
template <int dim>
void build_eigen_iterations
(std_cxx11::shared_ptr<EigenBase<dim> > eig_ptr, const ParameterHandler &prm)
{
  // TODO: we only have power iteration now, change later once we need to choose
  // different in group solvers
  eig_ptr =
  std_cxx11::shared_ptr<EigenBase<dim> >
  (new PowerIteration<dim> (prm));
}

/** \brief Function used to build pointer to instance of MGBase's derived class
 */
template <int dim>
void build_mg_iterations
(std_cxx11::shared_ptr<MGBase<dim> > mg_ptr, const ParameterHandler &prm)
{
  // TODO: fill this up once we have derived class of MGBase
  mg_ptr =
  std_cxx11::shared_ptr<MGBase<dim> > (new GaussSidel<dim> (prm));
}

/** \brief Function used to build pointer to instance of InGroupBase's derived class
 */
template <int dim>
void build_ig_iterations
(std_cxx11::shared_ptr<IGBase<dim> > ig_ptr, const ParameterHandler &prm)
{
  // TODO: we only have source iteration now, change later once we need to choose
  // different in group solvers
  bool do_nda = prm.get_bool ("do NDA");
  if (!do_nda)
    ig_ptr =
    std_cxx11::shared_ptr<IGBase<dim> > (new SourceIteration<dim> (prm));
}

/** \brief ConditionalOStream object used to output things on screen with processor 0
 */
ConditionalOStream pout
(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0);

#endif //__bart_tools_h__

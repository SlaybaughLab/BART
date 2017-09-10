#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>

#include "bart_tools.h"
#include ""
#include "../equations/even_parity.h"
#include "../aqdata/aq_lsgc"

template <int dim>
FE_Poly<TensorProductPolynomials<dim>,dim,dim>* build_finite_element
(ParameterHandler &prm)
{
  std::string discretization = prm.get ("spatial discretization");
  AssertThrow (discretization!="none",
               ExcMessage("valid spatial discretizations are dfem and cfem"))
  if (discretization=="dfem")
    return new FE_DGQ<dim> (p_order);
  else
    return new FE_Q<dim> (p_order);
}

std_cxx11::shared_ptr<PreconditionerSolver> build_linalg
(ParameterHandler &prm,
 std::string equation_name,
 unsigned int& n_total_vars)
{
  return std_cxx11::shared_ptr<PreconditionerSolver>
  (new PreconditionerSolver(prm, equation_name, n_total_vars));
}

template <int dim>
std_cxx11::shared_ptr<Iterations<dim> > build_iterations
(const ParameterHandler &prm,
 const DoFHandler<dim> &dof_handler,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
{
  return std_cxx11::shared_ptr<Iterations<dim> >
  (new Iterations<dim> (prm, msh_ptr, aqd_ptr, mat_ptr));
}

template <int dim>
std_cxx11::shared_ptr<MeshGenerator<dim> > build_mesh (ParameterHandler &prm)
{
  return std_cxx11::shared_ptr<MeshGenerator<dim> >
  (new MeshGenerator<dim>(prm));
}

std_cxx11::shared_ptr<MaterialProperties> build_material (ParameterHandler &prm)
{
  return std_cxx11::shared_ptr<MaterialProperties>
  (new MaterialProperties (prm));
}

template <int dim>
std_cxx11::shared_ptr<EquationBase<dim> > build_equation
(std::string space_angle_solver_name,
 ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
{
  // TODO: add NDA to it after having NDA class
  std_cxx11::shared_ptr<EquationBase<dim> > equation_pointer;
  if (transport_model_name=="ep")
    equation_pointer =
    std_cxx11::shared_ptr<EquationBase<dim> >
    (new EvenParity<dim> (prm, msh_ptr, aqd_ptr, mat_ptr));
  return equation_pointer;
}

template <int dim>
std_cxx11::shared_ptr<AQBase<dim> >
build_aq_model (ParameterHandler &prm)
{
  std::string aq_name = prm.get ("angular quadrature name");
  AssertThrow (aq_name!="none",
               ExcMessage("angular quadrature name incorrect or missing"));
  std_cxx11::shared_ptr<AQBase<dim> > aq_pointer;
  if (aq_name=="lsgc")
    aq_pointer = std_cxx11::shared_ptr<AQBase<dim> > (new LSGC<dim>(prm));
  return aq_pointer;
}

/** \brief Function used to build pointer to instance of InGroupBase's derived class
 */
template <int dim>
std_cxx11::shared_ptr<EigenBase<dim> > build_eigen_iterations (ParameterHandler &prm)
{
  // TODO: we only have power iteration now, change later once we need to choose
  // different in group solvers
  bool do_nda = prm.get_bool ("do NDA");
  if (!do_nda)
    return std_cxx11::shared_ptr<EigenBase<dim> >
    (new PowerIteration<dim> (prm));
}

/** \brief Function used to build pointer to instance of MGBase's derived class
 */
template <int dim>
std_cxx11::shared_ptr<MGBase<dim> > build_mg_iterations (ParameterHandler &prm)
{
  // TODO: fill this up once we have derived class of MGBase
}

/** \brief Function used to build pointer to instance of InGroupBase's derived class
 */
template <int dim>
std_cxx11::shared_ptr<IGBase<dim> > build_ig_iterations (ParameterHandler &prm)
{
  // TODO: we only have source iteration now, change later once we need to choose
  // different in group solvers
  bool do_nda = prm.get_bool ("do NDA");
  if (!do_nda)
    return std_cxx11::shared_ptr<IGBase<dim> >
    (new SourceIteration<dim> (prm));
}

void radio (std::string str)
{
  pout << str << std::endl;
}

void radio (std::string str1, std::string str2)
{
  pout << str1 << ": " << str2 << std::endl;
}

void radio (std::string str,
            double num)
{
  pout << str << ": " << num << std::endl;
}

void radio (std::string str1, unsigned int num1,
            std::string str2, unsigned int num2,
            std::string str3, unsigned int num3)
{
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
  {
    std::cout << str1 << ": " << num1 << ", ";
    std::cout << str2 << ": " << num2 << ", ";
    std::cout << str3 << ": " << num3 << std::endl;
  }
}

void radio (std::string str, unsigned int num)
{
  pout << str << ": " << num << std::endl;
}

void radio (std::string str, bool boolean)
{
  pout << str << ": " << (boolean?"true":"false") << std::endl;
}

void radio ()
{
  pout << "-------------------------------------" << std::endl << std::endl;
}

template std_cxx11::shared_ptr<Iterations<2> > build_iterations;
template std_cxx11::shared_ptr<Iterations<3> > build_iterations;
template std_cxx11::shared_ptr<EquationBase<2> > build_equation;
template std_cxx11::shared_ptr<EquationBase<3> > build_equation;
template std_cxx11::shared_ptr<AQBase<2> > build_aq_model;
template std_cxx11::shared_ptr<AQBase<3> > build_aq_model;
template std_cxx11::shared_ptr<MeshGenerator<2> > build_mesh;
template std_cxx11::shared_ptr<MeshGenerator<3> > build_mesh;
template std_cxx11::shared_ptr<EigenBase<2> > build_eigen_iterations;
template std_cxx11::shared_ptr<EigenBase<3> > build_eigen_iterations;
template std_cxx11::shared_ptr<MGBase<2> > build_mg_iterations;
template std_cxx11::shared_ptr<MGBase<3> > build_mg_iterations;
template std_cxx11::shared_ptr<InGroupBase<2> > build_ig_iterations;
template std_cxx11::shared_ptr<InGroupBase<3> > build_ig_iterations;


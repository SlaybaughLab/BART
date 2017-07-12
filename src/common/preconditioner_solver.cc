#include "../../include/common/preconditioner_solver.h"

PreconditionerSolver::PreconditionerSolver (ParameterHandler &prm,
                                            unsigned int& n_total_ho_vars,
                                            MPI_Comm &mpi_communicator)
:
n_group(prm.get_integer("number of groups")),
n_total_ho_vars(n_total_ho_vars),
mpi_communicator(mpi_communicator),
ho_linear_solver_name(prm.get("HO linear solver name")),
ho_preconditioner_name(prm.get("HO preconditioner name")),
do_nda(prm.get_bool("do NDA"))
{
	if (do_nda)
	{
		nda_linear_solver_name = prm.get ("NDA linear solver name");
		AssertThrow (nda_linear_solver_name!="cg",
			           ExcMessage("CG is prohibited for solving NDA equations"));
		nda_preconditioner_name = prm.get ("NDA preconditioner name");
	}
}

PreconditionerSolver::~PreconditionerSolver ()
{
}

// the following section is for HO solving/preconditioning
void PreconditionerSolver::initialize_ho_preconditioners
(std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
 std::vector<PETScWrappers::MPI::vector*> &ho_rhses)
{
  AssertThrow (n_total_ho_vars==ho_syses.size(),
               ExcMessage("num of HO system matrices should be equal to total variable number"));
  AssertThrow (n_total_ho_vars==ho_rhses.size(),
               ExcMessage("num of HO system rhs should be equal to total variable number"));
  if (ho_linear_solver_name!="direct")
  {
    linear_iters.resize (n_total_ho_vars);
    if (ho_preconditioner_name=="amg")
    {
      pre_ho_amg.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_amg[i] = (std_cxx11::shared_ptr<LA::MPI::PreconditionAMG>
                         (new LA::MPI::PreconditionAMG));
        LA::MPI::PreconditionAMG::AdditionalData data;
        if (transport_model_name=="fo" ||
            (transport_model_name=="ep" && have_reflective_bc))
          data.symmetric_operator = false;
        else
          data.symmetric_operator = true;
        pre_ho_amg[i]->initialize(*(ho_syses)[i], data);
      }
    }
    else if (ho_preconditioner_name=="bjacobi")
    {
      pre_ho_bjacobi.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_bjacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi>
        (new PETScWrappers::PreconditionBlockJacobi);
        pre_ho_bjacobi[i]->initialize(*(ho_syses)[i]);
      }
    }
    else if (ho_preconditioner_name=="jacobi")
    {
      pre_ho_jacobi.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_jacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi>
        (new PETScWrappers::PreconditionJacobi);
        pre_ho_jacobi[i]->initialize(*(ho_syses)[i]);
      }
    }
    else if (ho_preconditioner_name=="bssor")
    {
      pre_ho_eisenstat.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_eisenstat[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat>
        (new PETScWrappers::PreconditionEisenstat);
        PETScWrappers::PreconditionEisenstat::AdditionalData data(ssor_omega);
        pre_ho_eisenstat[i]->initialize(*(ho_syses)[i], data);
      }
    }
    else if (ho_preconditioner_name=="parasails")
    {
      pre_ho_parasails.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_parasails[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails>
                               (new PETScWrappers::PreconditionParaSails));
        if (transport_model_name=="fo" ||
            (transport_model_name=="ep" && have_reflective_bc))
        {
          PETScWrappers::PreconditionParaSails::AdditionalData data (2);
          pre_ho_parasails[i]->initialize(*(ho_syses)[i], data);
        }
        else
        {
          PETScWrappers::PreconditionParaSails::AdditionalData data (1);
          pre_ho_parasails[i]->initialize(*(ho_syses)[i], data);
        }
      }
    }
  }// not direct solver
  else
  {
    ho_direct.resize (n_total_ho_vars);
    direct_init = std::vector<bool> (n_total_ho_vars, false);
  }
  // initialize HO solver controls
  ho_cn.resize (n_total_ho_vars);
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
    ho_cn[i] = std_cxx11::shared_ptr<SolverControl>
    (new SolverControl(ho_rhses[i].size(),
                       1.0e-12*ho_rhses[i]->l1_norm()));
}

void PreconditionerSolver::ho_solve
(std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
 std::vector<PETScWrappers::MPI::vector*> &ho_rhses)
{
  AssertThrow (n_total_ho_vars==ho_syses.size(),
               ExcMessage("num of HO system matrices should be equal to total variable number"));
  AssertThrow (n_total_ho_vars==ho_rhses.size(),
               ExcMessage("num of HO system rhs should be equal to total variable number"));
  
}

// the following section is for NDA solving/preconditioning
void PreconditionerSolver::nda_solve
(std::vector<PETScWrappers::MPI::SparseMatrix*> &nda_syses,
 std::vector<PETScWrappers::MPI::Vector*> &nda_phis,
 std::vector<PETScWrappers::MPI::Vector*> &nda_rhses,
 unsigned int &g)
{
  // solve with specific solvers
  if (nda_linear_solver_name=="direct")
    nda_solve_direct (*nda_syses[g], *nda_phis, *nda_rhses, g);
  else if (nda_linear_solver_name=="bicgstab")
    nda_solve_bicgstab (*nda_syses, *nda_phis, *nda_rhses[g], g);
  else if (nda_linear_solver_name="gmres")
    nda_solve_gmres (*nda_syses[g], *nda_phis[g], *nda_rhses[g], g);
}

void PreconditionerSolver::reinit_nda_preconditioners
(std::vector<PETScWrappers::MPI::SparseMatrix*> &nda_syses)
{
  AssertThrow (nda_syses.size()==n_group,
               ExcMessage("There should be n_group NDA matrices"));
  if (nda_linear_solver_name=="direct")
  {
  }
  else
  {
  }
}

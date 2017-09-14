#include "preconditioner_solver.h"

PreconditionerSolver::PreconditionerSolver (const ParameterHandler &prm,
                                            std::string equation_name,
                                            unsigned int& n_total_vars)
:
equation_name(equation_name),
n_total_vars(n_total_vars)
{
  if (equation_name=="nda")
  {
    linear_solver_name = prm.get ("NDA linear solver name");
    AssertThrow (linear_solver_name!="cg",
                 ExcMessage("CG is prohibited for solving NDA equations"));
    preconditioner_name = prm.get ("NDA preconditioner name");
    if (preconditioner_name=="ssor")
      ssor_omega = prm.get_double ("NDA ssor factor");
  }
  else
  {
    linear_solver_name = prm.get ("HO linear solver name");
    if (equation_name=="ep")
    {
      have_reflective_bc = prm.get_bool ("have reflective BC");
      if (have_reflective_bc)
        AssertThrow (linear_solver_name!="cg",
                     ExcMessage("CG is prohibited for solving EP have reflective BC"));
      preconditioner_name = prm.get ("HO linear solver name");
      if (preconditioner_name=="ssor")
        ssor_omega = prm.get_double ("HO ssor factor");
    }
  }
}

PreconditionerSolver::~PreconditionerSolver ()
{
}

// the following section is for HO solving/preconditioning
void PreconditionerSolver::initialize_preconditioners
(std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
 std::vector<PETScWrappers::MPI::Vector*> &sys_rhses)
{
  AssertThrow (n_total_vars==sys_mats.size(),
               ExcMessage("num of system matrices should be equal to total variable number"));
  AssertThrow (n_total_vars==sys_rhses.size(),
               ExcMessage("num of system rhs should be equal to total variable number"));
  if (linear_solver_name!="direct")
  {
    linear_iters.resize (n_total_vars);
    if (preconditioner_name=="amg")
    {
      pre_amg.resize (n_total_vars);
      for (unsigned int i=0; i<n_total_vars; ++i)
      {
        pre_amg[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG>
                         (new PETScWrappers::PreconditionBoomerAMG));
        PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        if (equation_name=="fo" || equation_name=="nda" ||
            (equation_name=="ep" && have_reflective_bc))
          data.symmetric_operator = false;
        else
          data.symmetric_operator = true;
        pre_amg[i]->initialize(*sys_mats[i], data);
      }
    }
    else if (preconditioner_name=="bjacobi")
    {
      pre_bjacobi.resize (n_total_vars);
      for (unsigned int i=0; i<n_total_vars; ++i)
      {
        pre_bjacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi>
        (new PETScWrappers::PreconditionBlockJacobi);
        pre_bjacobi[i]->initialize(*sys_mats[i]);
      }
    }
    else if (preconditioner_name=="jacobi")
    {
      pre_jacobi.resize (n_total_vars);
      for (unsigned int i=0; i<n_total_vars; ++i)
      {
        pre_jacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi>
        (new PETScWrappers::PreconditionJacobi);
        pre_jacobi[i]->initialize(*sys_mats[i]);
      }
    }
    else if (preconditioner_name=="bssor")
    {
      pre_eisenstat.resize (n_total_vars);
      for (unsigned int i=0; i<n_total_vars; ++i)
      {
        pre_eisenstat[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat>
        (new PETScWrappers::PreconditionEisenstat);
        PETScWrappers::PreconditionEisenstat::AdditionalData data(ssor_omega);
        pre_eisenstat[i]->initialize(*sys_mats[i], data);
      }
    }
    else if (preconditioner_name=="parasails")
    {
      pre_parasails.resize (n_total_vars);
      for (unsigned int i=0; i<n_total_vars; ++i)
      {
        pre_parasails[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails>
                               (new PETScWrappers::PreconditionParaSails));
        if (equation_name=="fo" ||
            (equation_name=="ep" && have_reflective_bc))
        {
          PETScWrappers::PreconditionParaSails::AdditionalData data (2);
          pre_parasails[i]->initialize(*sys_mats[i], data);
        }
        else
        {
          PETScWrappers::PreconditionParaSails::AdditionalData data (1);
          pre_parasails[i]->initialize(*sys_mats[i], data);
        }
      }
    }
  }// not direct solver
  else
  {
    direct.resize (n_total_vars);
    direct_init = std::vector<bool> (n_total_vars, false);
  }
  // initialize HO solver controls
  cn.resize (n_total_vars);
  for (unsigned int i=0; i<n_total_vars; ++i)
    cn[i] = std_cxx11::shared_ptr<SolverControl>
    (new SolverControl(sys_rhses[i]->size(), 1.0e-12*sys_rhses[i]->l1_norm()));
}

void PreconditionerSolver::linear_algebra_solve
(std::vector<PETScWrappers::MPI::SparseMatrix*> &sys_mats,
 std::vector<PETScWrappers::MPI::Vector*> &sys_flxes,
 std::vector<PETScWrappers::MPI::Vector*> &sys_rhses,
 unsigned int &i/*component number*/)
{
  if (linear_solver_name=="cg")
  {
    PETScWrappers::SolverCG
    solver (*cn[i], MPI_COMM_WORLD);
    if (preconditioner_name=="amg")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_amg[i]);
    else if (preconditioner_name=="jacobi")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_jacobi[i]);
    else if (preconditioner_name=="bssor")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_eisenstat[i]);
    else if (preconditioner_name=="parasails")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_parasails[i]);
  }
  else if (linear_solver_name=="bicgstab")
  {
    PETScWrappers::SolverBicgstab
    solver (*cn[i], MPI_COMM_WORLD);
    if (preconditioner_name=="amg")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_amg[i]);
    else if (preconditioner_name=="jacobi")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_jacobi[i]);
    else if (preconditioner_name=="bssor")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_eisenstat[i]);
    else if (preconditioner_name=="parasails")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_parasails[i]);
  }
  else if (linear_solver_name=="gmres")
  {
    PETScWrappers::SolverGMRES
    solver (*cn[i], MPI_COMM_WORLD);
    if (preconditioner_name=="amg")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_amg[i]);
    else if (preconditioner_name=="jacobi")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_jacobi[i]);
    else if (preconditioner_name=="bssor")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_eisenstat[i]);
    else if (preconditioner_name=="parasails")
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *pre_parasails[i]);
  }
  else// if (linear_solver_name=="direct")
  {
    if (!direct_init[i])
    {
      direct[i] = std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS>
      (new PETScWrappers::SparseDirectMUMPS(*cn[i], MPI_COMM_WORLD));
      if (equation_name=="fo" ||
          (equation_name=="ep" && have_reflective_bc))
        direct[i]->set_symmetric_mode (false);
      else
        direct[i]->set_symmetric_mode (true);
      direct_init[i] = true;
    }
    direct[i]->solve (*sys_mats[i],
                      *sys_flxes[i],
                      *sys_rhses[i]);
  }
  // the linear_iters are for reporting linear solver status, test purpose only
  if (linear_solver_name!="direct")
    linear_iters[i] = cn[i]->last_step ();
}

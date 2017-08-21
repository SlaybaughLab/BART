#include "preconditioner_solver.h"

PreconditionerSolver::PreconditionerSolver (ParameterHandler &prm,
                                            unsigned int& n_total_ho_vars)
:
n_group(prm.get_integer("number of groups")),
n_total_ho_vars(n_total_ho_vars),
transport_model_name(prm.get("transport model")),
ho_linear_solver_name(prm.get("HO linear solver name")),
ho_preconditioner_name(prm.get("HO preconditioner name")),
do_nda(prm.get_bool("do NDA"))
{
  if (transport_model_name=="ep")
    have_reflective_bc = prm.get_bool ("have reflective BC");
  
  if (ho_preconditioner_name=="bssor")
    ho_ssor_omega = prm.get_double ("HO ssor factor");
  
	if (do_nda)
	{
		nda_linear_solver_name = prm.get ("NDA linear solver name");
		AssertThrow (nda_linear_solver_name!="cg",
			           ExcMessage("CG is prohibited for solving NDA equations"));
		nda_preconditioner_name = prm.get ("NDA preconditioner name");
    if (nda_preconditioner_name=="ssor")
      nda_ssor_omega = prm.get_double ("NDA ssor factor");
	}
}

PreconditionerSolver::~PreconditionerSolver ()
{
}

// the following section is for HO solving/preconditioning
void PreconditionerSolver::initialize_ho_preconditioners
(std::vector<PETScWrappers::MPI::SparseMatrix*> &ho_syses,
 std::vector<PETScWrappers::MPI::Vector*> &ho_rhses)
{
  AssertThrow (n_total_ho_vars==ho_syses.size(),
               ExcMessage("num of HO system matrices should be equal to total variable number"));
  AssertThrow (n_total_ho_vars==ho_rhses.size(),
               ExcMessage("num of HO system rhs should be equal to total variable number"));
  if (ho_linear_solver_name!="direct")
  {
    ho_linear_iters.resize (n_total_ho_vars);
    if (ho_preconditioner_name=="amg")
    {
      pre_ho_amg.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_amg[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG>
                         (new PETScWrappers::PreconditionBoomerAMG));
        PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        if (transport_model_name=="fo" ||
            (transport_model_name=="ep" && have_reflective_bc))
          data.symmetric_operator = false;
        else
          data.symmetric_operator = true;
        pre_ho_amg[i]->initialize(*ho_syses[i], data);
      }
    }
    else if (ho_preconditioner_name=="bjacobi")
    {
      pre_ho_bjacobi.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_bjacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi>
        (new PETScWrappers::PreconditionBlockJacobi);
        pre_ho_bjacobi[i]->initialize(*ho_syses[i]);
      }
    }
    else if (ho_preconditioner_name=="jacobi")
    {
      pre_ho_jacobi.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_jacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi>
        (new PETScWrappers::PreconditionJacobi);
        pre_ho_jacobi[i]->initialize(*ho_syses[i]);
      }
    }
    else if (ho_preconditioner_name=="bssor")
    {
      pre_ho_eisenstat.resize (n_total_ho_vars);
      for (unsigned int i=0; i<n_total_ho_vars; ++i)
      {
        pre_ho_eisenstat[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat>
        (new PETScWrappers::PreconditionEisenstat);
        PETScWrappers::PreconditionEisenstat::AdditionalData data(ho_ssor_omega);
        pre_ho_eisenstat[i]->initialize(*ho_syses[i], data);
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
          pre_ho_parasails[i]->initialize(*ho_syses[i], data);
        }
        else
        {
          PETScWrappers::PreconditionParaSails::AdditionalData data (1);
          pre_ho_parasails[i]->initialize(*ho_syses[i], data);
        }
      }
    }
  }// not direct solver
  else
  {
    ho_direct.resize (n_total_ho_vars);
    ho_direct_init = std::vector<bool> (n_total_ho_vars, false);
  }
  // initialize HO solver controls
  ho_cn.resize (n_total_ho_vars);
  for (unsigned int i=0; i<n_total_ho_vars; ++i)
    ho_cn[i] = std_cxx11::shared_ptr<SolverControl>
    (new SolverControl(ho_rhses[i]->size(), 1.0e-12*ho_rhses[i]->l1_norm()));
}

void PreconditionerSolver::ho_solve
(PETScWrappers::MPI::SparseMatrix &ho_sys,
 PETScWrappers::MPI::Vector &ho_psi,
 PETScWrappers::MPI::Vector &ho_rhs,
 unsigned int i/*component number*/)
{
  if (ho_linear_solver_name=="cg")
  {
    PETScWrappers::SolverCG
    solver (*ho_cn[i], MPI_COMM_WORLD);
    if (ho_preconditioner_name=="amg")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_amg[i]);
    else if (ho_preconditioner_name=="jacobi")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_jacobi[i]);
    else if (ho_preconditioner_name=="bssor")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_eisenstat[i]);
    else if (ho_preconditioner_name=="parasails")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_parasails[i]);
  }
  else if (ho_linear_solver_name=="bicgstab")
  {
    PETScWrappers::SolverBicgstab
    solver (*ho_cn[i], MPI_COMM_WORLD);
    if (ho_preconditioner_name=="amg")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_amg[i]);
    else if (ho_preconditioner_name=="jacobi")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_jacobi[i]);
    else if (ho_preconditioner_name=="bssor")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_eisenstat[i]);
    else if (ho_preconditioner_name=="parasails")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_parasails[i]);
  }
  else if (ho_linear_solver_name=="gmres")
  {
    PETScWrappers::SolverGMRES
    solver (*ho_cn[i], MPI_COMM_WORLD);
    if (ho_preconditioner_name=="amg")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_amg[i]);
    else if (ho_preconditioner_name=="jacobi")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_jacobi[i]);
    else if (ho_preconditioner_name=="bssor")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_eisenstat[i]);
    else if (ho_preconditioner_name=="parasails")
      solver.solve (ho_sys,
                    ho_psi,
                    ho_rhs,
                    *pre_ho_parasails[i]);
  }
  else// if (linear_solver_name=="direct")
  {
    if (!ho_direct_init[i])
    {
      ho_direct[i] = std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS>
      (new PETScWrappers::SparseDirectMUMPS(*ho_cn[i], MPI_COMM_WORLD));
      if (transport_model_name=="fo" ||
          (transport_model_name=="ep" && have_reflective_bc))
        ho_direct[i]->set_symmetric_mode (false);
      else
        ho_direct[i]->set_symmetric_mode (true);
      ho_direct_init[i] = true;
    }
    ho_direct[i]->solve (ho_sys,
                         ho_psi,
                         ho_rhs);
  }
  // the ho_linear_iters are for reporting linear solver status, test purpose only
  if (ho_linear_solver_name!="direct")
    ho_linear_iters[i] = ho_cn[i]->last_step ();
}

// the following section is for NDA solving/preconditioning
// Unlike HO system, preconditioner will be reinit every outer
// iteration
void PreconditionerSolver::reinit_nda_preconditioners
(std::vector<PETScWrappers::MPI::SparseMatrix*> &nda_syses,
 std::vector<PETScWrappers::MPI::Vector*> &nda_rhses)
{
  AssertThrow (nda_syses.size()==n_group,
               ExcMessage("There sndauld be n_group NDA matrices"));
  AssertThrow (nda_rhses.size()==n_group,
               ExcMessage("num of NDA system rhs sndauld be equal to total group num"));
  if (nda_linear_solver_name!="direct")
  {
    nda_linear_iters.resize (n_group);
    if (nda_preconditioner_name=="amg")
    {
      pre_nda_amg.resize (n_group);
      for (unsigned int i=0; i<n_group; ++i)
      {
        pre_nda_amg[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionBoomerAMG>
                          (new PETScWrappers::PreconditionBoomerAMG));
        PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        data.symmetric_operator = false;
        pre_nda_amg[i]->initialize(*nda_syses[i], data);
      }
    }
    else if (nda_preconditioner_name=="bjacobi")
    {
      pre_nda_bjacobi.resize (n_group);
      for (unsigned int i=0; i<n_group; ++i)
      {
        pre_nda_bjacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionBlockJacobi>
        (new PETScWrappers::PreconditionBlockJacobi);
        pre_nda_bjacobi[i]->initialize(*nda_syses[i]);
      }
    }
    else if (nda_preconditioner_name=="jacobi")
    {
      pre_nda_jacobi.resize (n_group);
      for (unsigned int i=0; i<n_group; ++i)
      {
        pre_nda_jacobi[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionJacobi>
        (new PETScWrappers::PreconditionJacobi);
        pre_nda_jacobi[i]->initialize(*nda_syses[i]);
      }
    }
    else if (nda_preconditioner_name=="bssor")
    {
      pre_nda_eisenstat.resize (n_group);
      for (unsigned int i=0; i<n_group; ++i)
      {
        pre_nda_eisenstat[i] = std_cxx11::shared_ptr<PETScWrappers::PreconditionEisenstat>
        (new PETScWrappers::PreconditionEisenstat);
        PETScWrappers::PreconditionEisenstat::AdditionalData data(nda_ssor_omega);
        pre_nda_eisenstat[i]->initialize(*nda_syses[i], data);
      }
    }
    else if (nda_preconditioner_name=="parasails")
    {
      pre_nda_parasails.resize (n_group);
      for (unsigned int i=0; i<n_group; ++i)
      {
        pre_nda_parasails[i] = (std_cxx11::shared_ptr<PETScWrappers::PreconditionParaSails>
                                (new PETScWrappers::PreconditionParaSails));
        // set the symmetric pattern to false
        PETScWrappers::PreconditionParaSails::AdditionalData data (2);
        pre_nda_parasails[i]->initialize(*nda_syses[i], data);
      }
    }
  }// not direct solver
  else
  {
    nda_direct.resize (n_group);
    nda_direct_init = std::vector<bool> (n_group, false);
  }
  // initialize nda solver controls
  nda_cn.resize (n_group);
  for (unsigned int i=0; i<n_group; ++i)
    nda_cn[i] = std_cxx11::shared_ptr<SolverControl>
    (new SolverControl(nda_rhses[i]->size(),
                       1.0e-12*nda_rhses[i]->l1_norm()));
}

void PreconditionerSolver::nda_solve
(PETScWrappers::MPI::SparseMatrix &nda_sys,
 PETScWrappers::MPI::Vector &nda_phi,
 PETScWrappers::MPI::Vector &nda_rhs,
 unsigned int &g)
{
  AssertThrow (nda_linear_solver_name!="cg",
               ExcMessage("CG is not available to NDA solving"));
  // solve with specific solvers
  if (nda_linear_solver_name=="direct")
  {
    if (!nda_direct_init[g])
    {
      nda_direct[g] = std_cxx11::shared_ptr<PETScWrappers::SparseDirectMUMPS>
      (new PETScWrappers::SparseDirectMUMPS(*nda_cn[g], MPI_COMM_WORLD));
      nda_direct[g]->set_symmetric_mode (false);
      nda_direct_init[g] = true;
    }
    nda_direct[g]->solve (nda_sys,
                          nda_phi,
                          nda_rhs);
  }
  else if (nda_linear_solver_name=="bicgstab")
  {
    PETScWrappers::SolverBicgstab
    solver (*nda_cn[g], MPI_COMM_WORLD);
    if (nda_preconditioner_name=="amg")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_amg[g]);
    else if (nda_preconditioner_name=="parasails")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_parasails[g]);
    else if (nda_preconditioner_name=="bssor")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_eisenstat[g]);
    else if (nda_preconditioner_name=="bjacobi")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_bjacobi[g]);
    else if (nda_preconditioner_name=="jacobi")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_jacobi[g]);
  }
  else if (nda_linear_solver_name=="gmres")
  {
    PETScWrappers::SolverGMRES
    solver (*nda_cn[g], MPI_COMM_WORLD);
    if (nda_preconditioner_name=="amg")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_amg[g]);
    else if (nda_preconditioner_name=="parasails")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_parasails[g]);
    else if (nda_preconditioner_name=="bssor")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_eisenstat[g]);
    else if (nda_preconditioner_name=="bjacobi")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_bjacobi[g]);
    else if (nda_preconditioner_name=="jacobi")
      solver.solve (nda_sys,
                    nda_phi,
                    nda_rhs,
                    *pre_nda_jacobi[g]);
  }
}


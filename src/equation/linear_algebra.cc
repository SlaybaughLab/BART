#include "linear_algebra.h"

LinearAlgebra::LinearAlgebra(const dealii::ParameterHandler &prm,
    const std::string &equ_name,
    const int &n_total_vars)
    :
    equ_name_(equ_name),
    n_total_vars_(n_total_vars),
    have_reflective_bc_(prm.get_bool("have reflective boundary")) {
  if (equ_name_=="nda" || equ_name_=="tg_nda") {
    linear_solver_name_ = prm.get ("nda linear solver name");
    preconditioner_name_ = prm.get ("nda preconditioner name");
    if (preconditioner_name_=="ssor")
      ssor_omega_ = prm.get_double ("nda ssor factor");
  } else {
    linear_solver_name_ = prm.get ("ho linear solver name");
    preconditioner_name_ = prm.get ("ho linear solver name");
    if (preconditioner_name_=="ssor")
      ssor_omega_ = prm.get_double ("ho ssor factor");
  }
}

LinearAlgebra::~LinearAlgebra () {}

void LinearAlgebra::InitPrecond (
    std::unordered_map<int, dealii::PETScWrappers::MPI::SparseMatrix*> &sys_mats,
    std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> &sys_rhses) {
  AssertThrow (n_total_vars_ == static_cast<int>(sys_mats.size()),
      dealii::ExcMessage("num of system matrices should be equal to total variable number"));
  if (linear_solver_name_!="direct") {
    if (preconditioner_name_=="amg") {
      for (int i=0; i<n_total_vars_; ++i) {
        p_amg_[i] =
            std::unique_ptr<dealii::PETScWrappers::PreconditionBoomerAMG> (
            new dealii::PETScWrappers::PreconditionBoomerAMG);
        dealii::PETScWrappers::PreconditionBoomerAMG::AdditionalData data;
        if (equ_name_=="fo" || equ_name_=="nda" ||
            (equ_name_=="ep" && have_reflective_bc_))
          data.symmetric_operator = false;
        else
          data.symmetric_operator = true;
        p_amg_[i]->initialize(*sys_mats[i], data);
      }
    } else if (preconditioner_name_=="bjacobi") {
      for (int i=0; i<n_total_vars_; ++i) {
        p_bjacobi_[i] =
          std::unique_ptr<dealii::PETScWrappers::PreconditionBlockJacobi> (
          new dealii::PETScWrappers::PreconditionBlockJacobi);
        p_bjacobi_[i]->initialize(*sys_mats[i]);
      }
    } else if (preconditioner_name_=="jacobi") {
      for (int i=0; i<n_total_vars_; ++i) {
        p_jacobi_[i] =
            std::unique_ptr<dealii::PETScWrappers::PreconditionJacobi> (
            new dealii::PETScWrappers::PreconditionJacobi);
        p_jacobi_[i]->initialize(*sys_mats[i]);
      }
    } else if (preconditioner_name_=="bssor") {
      for (int i=0; i<n_total_vars_; ++i) {
        p_eisenstat_[i] =
            std::unique_ptr<dealii::PETScWrappers::PreconditionEisenstat> (
            new dealii::PETScWrappers::PreconditionEisenstat);
        dealii::PETScWrappers::PreconditionEisenstat::AdditionalData data(ssor_omega_);
        p_eisenstat_[i]->initialize(*sys_mats[i], data);
      }
    } else if (preconditioner_name_=="parasails") {
      for (int i=0; i<n_total_vars_; ++i) {
        p_parasails_[i] =
            std::unique_ptr<dealii::PETScWrappers::PreconditionParaSails> (
            new dealii::PETScWrappers::PreconditionParaSails);
        if (equ_name_=="fo" || equ_name_=="nda" || equ_name_=="tg_nda" ||
            (equ_name_=="ep" && have_reflective_bc_)) {
          dealii::PETScWrappers::PreconditionParaSails::AdditionalData data (2);
          p_parasails_[i]->initialize(*sys_mats[i], data);
        } else {
          dealii::PETScWrappers::PreconditionParaSails::AdditionalData data (1);
          p_parasails_[i]->initialize(*sys_mats[i], data);
        }
      }
    }
  } else {
    for (int i=0; i<n_total_vars_; ++i)
      direct_init_[i] = false;
  }
  // TODO: find a better way for solver control
  // initialize solver control
  cn_ = std::unique_ptr<dealii::SolverControl> (
      new dealii::SolverControl(sys_rhses[0]->size(),
      sys_rhses[0]->size()*1.0e-18));
}

void LinearAlgebra::LinAlgSolve (
    std::unordered_map<int, dealii::PETScWrappers::MPI::SparseMatrix*> &sys_mats,
    std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> &sys_flxes,
    std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> &sys_rhses,
    std::unordered_map<int, dealii::ConstraintMatrix*> &,
    const int &i) {
  if (linear_solver_name_=="cg") {
    dealii::PETScWrappers::SolverCG solver (*cn_, MPI_COMM_WORLD);
    if (preconditioner_name_=="amg") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_amg_[i]);
    } else if (preconditioner_name_=="jacobi") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_jacobi_[i]);
    } else if (preconditioner_name_=="bssor") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_eisenstat_[i]);
    } else if (preconditioner_name_=="parasails") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_parasails_[i]);
    }
  } else if (linear_solver_name_=="bicgstab") {
    dealii::PETScWrappers::SolverBicgstab solver (*cn_, MPI_COMM_WORLD);
    if (preconditioner_name_=="amg") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_amg_[i]);
    } else if (preconditioner_name_=="jacobi") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_jacobi_[i]);
    } else if (preconditioner_name_=="bssor") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_eisenstat_[i]);
    } else if (preconditioner_name_=="parasails") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_parasails_[i]);
    }
  } else if (linear_solver_name_=="gmres") {
    dealii::PETScWrappers::SolverGMRES solver (*cn_, MPI_COMM_WORLD);
    if (preconditioner_name_=="amg") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_amg_[i]);
    } else if (preconditioner_name_=="jacobi") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_jacobi_[i]);
    } else if (preconditioner_name_=="bssor") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_eisenstat_[i]);
    } else if (preconditioner_name_=="parasails") {
      solver.solve (*sys_mats[i],
                    *sys_flxes[i],
                    *sys_rhses[i],
                    *p_parasails_[i]);
    }
  } else {
    // if the solvers have not been initilized yet
    if (!direct_init_[i]) {
      direct_[i] = std::unique_ptr<dealii::PETScWrappers::SparseDirectMUMPS> (
          new dealii::PETScWrappers::SparseDirectMUMPS(*cn_, MPI_COMM_WORLD));
      if (equ_name_=="fo" || equ_name_=="nda" || equ_name_=="tg_nda" ||
          (equ_name_=="ep" && have_reflective_bc_)) {
        direct_[i]->set_symmetric_mode (false);
      } else {
        direct_[i]->set_symmetric_mode (true);
      }
      direct_init_[i] = true;
    }
    direct_[i]->solve (*sys_mats[i], *sys_flxes[i], *sys_rhses[i]);
  }
  // TODO: apply constraints in future projects
  //constraints.distribute (*sys_flxes[i]);
  // the linear_iters are for reporting linear solver status, test purpose only
  if (linear_solver_name_!="direct")
    linear_iters_[i] = cn_->last_step ();
}

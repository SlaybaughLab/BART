#include "bart_driver.h"

template <int dim>
BARTDriver<dim>::BARTDriver (dealii::ParameterHandler &prm)
    :
    n_proc_(dealii::Utilities::MPI::n_processors),
    tria(typename dealii::Triangulation<dim>::MeshSmoothing
         dealii::Triangulation<dim>::smoothing_on_refinement |
         dealii::Triangulation<dim>::smoothing_on_coarsening)),
    distributed_tria(MPI_COMM_WORLD,
        typename dealii::Triangulation<dim>::MeshSmoothing
        dealii::Triangulation<dim>::smoothing_on_refinement |
        dealii::Triangulation<dim>::smoothing_on_coarsening)),
    is_eigen_problem_(prm.get_bool(prm.get_bool("do eigenvalue calculations"))),
    do_nda_(prm.get_bool("do NDA")),
    ho_equ_name_(prm.get("transport model")),
    dat_ptr_(std::shared_ptr<FundamentalData<dim>>(
        new FundamentalData<dim>(prm, n_proc_==1?tria:distributed_tria))),
    iter_ptr_(std::shared_ptr<Iterations<dim>>(new Iterations<dim>(prm))),
    equ_ptrs_(bbuilders::BuildEquations(prm, dat_ptr_)) {}

template <int dim>
BARTDriver<dim>::~BARTDriver () {}

template <int dim>
void BARTDriver<dim>::MakeGrid() {
  dat_ptr_->mesh->MakeGrid(n_proc_==1?tria:distributed_tria);
}

template <int dim>
void BARTDriver<dim>::InitMatVec() {
  //TODO: the following is assuming HO and LO are using the same finite elements
  //s.t. only one DoFHandler object is necessary. Fix this in future.
  dat_ptr_->dof_handler.distribute_dofs (*(dat_ptr_->fe_data->fe[ho_equ_name_]));
  local_owned_dofs_ = dat_ptr_->dof_handler.locally_owned_dofs ();
  dealii::DoFTools::extract_locally_relevant_dofs (dat_ptr_->dof_handler,
      local_relevant_dofs_);

  //TODO: in future if AMR is of interest, dummy_constraints_ has to be changed
  //to correct way.
  dummy_constraints_.clear ();
  dummy_constraints_.reinit (local_relevant_dofs_);
  dealii::DoFTools::make_hanging_node_constraints (dat_ptr_->dof_handler,
      dummy_constraints_);
  dummy_constraints_.close ();

  dealii::DynamicSparsityPattern dsp (local_relevant_dofs_);

  //TODO: current implementation requires HO and LO uses same discretization.
  //This has to be changed when using
  AssertThrow(ho_discretization_==nda_discretization_,
      dealii::ExcMessage("HO and NDA has to be using the same discretization for now"));
  if (ho_discretization=="dfem")
    dealii::DoFTools::make_flux_sparsity_pattern (dat_ptr_->dof_handler,
        dsp, dummy_constraints_, false);
  else
    dealii::DoFTools::make_sparsity_pattern (dat_ptr_->dof_handler,
        dsp, dummy_constraints_, false);

  // setting up dsp with telling communicator and relevant dofs
  dealii::SparsityTools::distribute_sparsity_pattern (dsp,
      dat_ptr_->dof_handler.n_locally_owned_dofs_per_processor (),
      MPI_COMM_WORLD, local_relevant_dofs_);

  const int n_total_vars = ho_equ_name_=="diffusion"?n_group_:
      dat_ptr_->aq->GetNTotalVars();

  std::shared_ptr<MatrixVector> mv = dat_ptr_->mat_vec;
  std::unordered_map<int, dealii::PETScWrappers::MPI::SparseMatrix*> mp_mat;
  std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> mp_vec;
  std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> mp_rhs;
  std::unordered_map<int, dealii::PETScWrappers::MPI::Vector*> mp_fix;
  for (int i=0; i<n_total_vars; ++i) {
    mp_mat[i] = new PETScWrappers::MPI::SparseMatrix;
    mp_mat[i]->reinit (local_owned_dofs_, local_owned_dofs_, dsp, MPI_COMM_WORLD);
    //TODO: relevant dofs are needed if adaptive mesh refinement is of interest
    mp_vec[i] = new PETScWrappers::MPI::Vector;
    mp_vec[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
    mp_rhs[i] = new PETScWrappers::MPI::Vector;
    mp_rhs[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
    mp_fix[i] = new PETScWrappers::MPI::Vector;
    mp_fix[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
  }
  mv->sys_mats[ho_equ_name_] = mp_mat;
  mp_mat.clear ();
  mv->sys_flxes[ho_equ_name_] = mp_vec;
  mp_vec.clear ();
  mv->sys_rhses[ho_equ_name_] = mp_rhs;
  mp_rhs.clear ();
  mv->sys_fixed_rhses[ho_equ_name_] = mp_fix;
  mp_fix.clear ();

  if (do_nda_) {
    n_total_vars = n_group_;
    for (int i=0; i<n_total_vars; ++i) {
      mp_mat[i] = new PETScWrappers::MPI::SparseMatrix;
      mp_mat[i]->reinit (local_owned_dofs_, local_owned_dofs_, dsp, MPI_COMM_WORLD);
      //TODO: relevant dofs are needed if adaptive mesh refinement is of interest
      mp_vec[i] = new PETScWrappers::MPI::Vector;
      mp_vec[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
      mp_rhs[i] = new PETScWrappers::MPI::Vector;
      mp_rhs[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
      mp_fix[i] = new PETScWrappers::MPI::Vector;
      mp_fix[i]->reinit (local_owned_dofs_, MPI_COMM_WORLD);
    }

    mv->sys_mats["nda"] = mp_mat;
    mp_mat.clear ();
    mv->sys_flxes["nda"] = mp_vec;
    mp_vec.clear ();
    mv->sys_rhses["nda"] = mp_rhs;
    mp_rhs.clear ();
    mv->sys_fixed_rhses["nda"] = mp_fix;
    mp_fix.clear ();

    //TODO: need to add one more equation if two-grid energy acceleration is
    //used for NDA
  }

  // TODO: fix the following part for Anisotropic scattering
  std::unordered_set<std::string> nm = do_nda_?{ho_equ_name_,"nda"}:
      {ho_equ_name_};
  for (const auto &name : nm) {
    std::map<std::tuple<int,int,int>, dealii::Vector<double>> mp;
    for (int g=0; g<n_group_; ++g) {
      mp[std::make_tuple(g,0,0)] = dealii::Vector(*(mv->sys_flxes[ho_equ_name_][0]));
      // set values to be 1.
      mp[std::make_tuple(g,0,0)] = 1.0;
    }
    mv->moments[name] = mp;
  }
}

template <int dim>
void BARTDriver<dim>::DoIterations() {
  iter_cls_.DoIterations(equ_ptrs_);
}

template <int dim>
void BARTDriver<dim>::OutputResults () const {
  dat_ptr_->pcout << "output results" << std::endl;
  dealii::DataOut<dim> data_out;
  data_out.attach_dof_handler (dat_ptr_->dof_handler);

  for (auto &mmts : dat_ptr_->moments) {
    for (int g=0; g<n_group_; ++g) {
      std::ostringstream os;
      os << mmts.first << "-phi-" << g;
      data_out.add_data_vector (mmts.second[std::make_tuple(g,0,0)], os.str ());
    }
  }

  dealii::Vector<float> subdomain (
      n_proc_==1?tria.n_active_cells ():distributed_tria.n_active_cells());
  const int proc_id = n_proc_==1?0:distributed_tria.locally_owned_subdomain ();
  for (int i=0; i<subdomain.size(); ++i)
    subdomain(i) = proc_id;
  data_out.add_data_vector (subdomain, "subdomain");

  if (is_eigen_problem_) {
    // add keff to the output file
    dealii::Vector<float> keffs (
        n_proc_==1?tria.n_active_cells ():distributed_tria.n_active_cells());
    const double keff = iter_ptr_->GetKeff();
    for (int i=0; i<keffs.size(); ++i)
      keffs(i) = keff;
    data_out.add_data_vector (keffs, "keff");
  }

  data_out.build_patches ();

  const std::string local_fname = output_fname_ +
      "-" + dealii::Utilities::int_to_string(proc_id, 4);
  std::ofstream output ((local_fname + ".vtu").c_str ());
  data_out.write_vtu (output);

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0) {
    std::vector<std::string> filenames;
    for (int i=0; i<n_proc_; ++i)
      filenames.push_back (local_fname);
    std::ostringstream os;
    os << output_fname_ << ".pvtu";
    std::ofstream master_output ((os.str()).c_str ());
    data_out.write_pvtu_record (master_output, filenames);
  }
}

template <int dim>
void BARTDriver<dim>::DriveBART () {
  // produce a grid
  MakeGrid();
  // initialize PETSc data structures
  InitMatVec();
  // Perform iterations
  DoIterations();
  // Output the results
  OutputResults ();
}

template class BARTDriver<1>;
template class BARTDriver<2>;
template class BARTDriver<3>;

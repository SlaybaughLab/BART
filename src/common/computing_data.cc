#include "computing_data.h"
#include "bart_builder.h"

#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>

XSections::XSections (MaterialBase& material)
    :
    sigt(material.GetSigT()),
    inv_sigt(material.GetInvSigT()),
    q(material.GetQ()),
    q_per_ster(material.GetQPerSter()),
    is_material_fissile(material.GetFissileIDMap()),
    nu_sigf(material.GetNuSigF()),
    sigs(material.GetSigS()),
    sigs_per_ster(material.GetSigSPerSter()),
    fiss_transfer(material.GetChiNuSigF()),
    fiss_transfer_per_ster(material.GetChiNuSigFPerSter())
{}

template <int dim>
FundamentalData<dim>::FundamentalData (dealii::ParameterHandler &prm,
    dealii::Triangulation<dim> &tria)
    :
    pcout(std::cout,
        (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)),
    aq(AQBase<dim>::CreateAQ(prm)),
    //material(prm),
    mesh(prm),
    mat_vec(std::shared_ptr<MatrixVector> (new MatrixVector())),
    //xsec(std::shared_ptr<XSections> (new XSections(material))),
    fe_data(prm),
    dof_handler(tria) {
  MaterialProtobuf material{prm};
  xsec = std::make_shared<XSections>(material);
  aq->MakeAQ();
}

template <int dim>
FundamentalData<dim>::~FundamentalData () {}

template <int dim>
FEData<dim>::FEData (const dealii::ParameterHandler &prm) {
  bbuilders::BuildFESpaces(prm, fe);
  const dealii::UpdateFlags update_flags =
      dealii::update_values | dealii::update_gradients |
      dealii::update_quadrature_points |
      dealii::update_JxW_values;
  const dealii::UpdateFlags update_face_flags = dealii::update_normal_vectors |
      dealii::update_values | dealii::update_gradients |
      dealii::update_quadrature_points |
      dealii::update_JxW_values;

  for (auto &f : fe) {
    // f.first is the equation name, f.second is the correspoding FiniteElement
    // object
    // TODO: all fe use the same finite element for now.
    if (f.first!="nda" && f.first!="tg_nda")
      discretization[f.first] = prm.get("ho spatial discretization");
    else
      discretization[f.first] = prm.get("nda spatial discretization");
    p_order[f.first] = prm.get_integer("finite element polynomial degree");
    q_rule[f.first] = std::shared_ptr<dealii::QGauss<dim>>(
        new dealii::QGauss<dim>(p_order[f.first]+1));
    qf_rule[f.first] = std::shared_ptr<dealii::QGauss<dim-1>>(
        new dealii::QGauss<dim-1>(p_order[f.first]+1));

    fv[f.first] = std::shared_ptr<dealii::FEValues<dim>> (
        new dealii::FEValues<dim> (*f.second, *q_rule[f.first], update_flags));

    fvf[f.first] = std::shared_ptr<dealii::FEFaceValues<dim>> (
        new dealii::FEFaceValues<dim>(
            *f.second, *qf_rule[f.first], update_face_flags));

    fvf_nei[f.first] = std::shared_ptr<dealii::FEFaceValues<dim>> (
        new dealii::FEFaceValues<dim>(
            *f.second, *qf_rule[f.first], update_face_flags));

    dofs_per_cell[f.first] = fe[f.first]->dofs_per_cell;
    n_q[f.first] = q_rule[f.first]->size();
    n_qf[f.first] = qf_rule[f.first]->size();

    local_dof_indices[f.first].resize (dofs_per_cell[f.first]);
    neigh_dof_indices[f.first].resize (dofs_per_cell[f.first]);

    if (f.first=="nda") {
      // corrections related objects
      q_rule["corr"] = std::shared_ptr<dealii::QGauss<dim>>(
          new dealii::QGauss<dim>(p_order[f.first]+3));
      qf_rule["corr"] = std::shared_ptr<dealii::QGauss<dim-1>>(
          new dealii::QGauss<dim-1>(p_order[f.first]+3));

      fv["corr"] = std::shared_ptr<dealii::FEValues<dim>> (
          new dealii::FEValues<dim> (*f.second, *q_rule["corr"], update_face_flags));

      fvf["corr"] = std::shared_ptr<dealii::FEFaceValues<dim>> (
          new dealii::FEFaceValues<dim>(
              *f.second, *qf_rule["corr"], update_face_flags));

      fvf_nei["corr"] = std::shared_ptr<dealii::FEFaceValues<dim>> (
          new dealii::FEFaceValues<dim>(
              *f.second, *qf_rule["corr"], update_face_flags));

      n_q["corr"] = q_rule["corr"]->size();
      n_qf["corr"] = qf_rule["corr"]->size();
    }
  }
}

template <int dim>
void FEData<dim>::InitializeFESpaces(const dealii::ParameterHandler &prm) {
    // getting parameter values
  const bool do_nda = prm.get_bool ("do nda");
  const int p_order = prm.get_integer("finite element polynomial degree");
  const std::string ho_discretization = prm.get ("ho spatial discretization");
  const std::string nda_discretization = prm.get ("nda spatial discretization");
  const std::string ho_equ_name = prm.get ("transport model");

  std::unordered_map<std::string, int> discretization_ind = {
      {"cfem", 0}, {"dfem", 1}, {"cmfd", 2}, {"rtk", 3}};

  switch (discretization_ind[ho_discretization]) {
    case 0:
      fe[ho_equ_name] = new dealii::FE_Q<dim> (p_order);
      break;

    case 1:
      fe[ho_equ_name] = new dealii::FE_DGQ<dim> (p_order);
      break;

    default:
      AssertThrow (false,
          dealii::ExcMessage("Invalid HO discretization name"));
      break;
  }

  if (do_nda) {
    switch (discretization_ind[nda_discretization]) {
      case 0:
        fe["nda"] = new dealii::FE_Q<dim> (p_order);
        break;

      case 1:
        fe["nda"] = new dealii::FE_DGQ<dim> (p_order);
        break;

      case 2:
        fe["nda"] = new dealii::FE_DGQ<dim> (0);
        break;

      case 3:
        fe["nda"] = new dealii::FE_RaviartThomas<dim> (p_order);
        break;

      default:
        AssertThrow (false,
            dealii::ExcMessage("Invalid NDA discretization name"));
        break;
    }
  }
}


template struct FundamentalData<1>;
template struct FundamentalData<2>;
template struct FundamentalData<3>;

template struct FEData<1>;
template struct FEData<2>;
template struct FEData<3>;

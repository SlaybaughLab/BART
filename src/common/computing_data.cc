#include "computing_data.h"
#include "bart_builder.h"

#include <deal.II/fe/fe_update_flags.h>

// XSections::XSections (Materials& material)
//     :
//     sigt(material.GetSigT()),
//     inv_sigt(material.GetInvSigT()),
//     q(material.GetQ()),
//     q_per_ster(material.GetQPerSter()),
//     is_material_fissile(material.GetFissileIDMap()),
//     nu_sigf(material.GetNuSigf()),
//     sigs(material.GetSigS()),
//     sigs_per_ster(material.GetSigSPerSter()),
//     fiss_transfer(material.GetFissTransfer()),
//     fiss_transfer_per_ster(material.GetFissTransferPerSter()) {}

XSections::XSections (MaterialPropertiesI& material_properties)
    :
    sigt(material_properties.GetSigT()),
    inv_sigt(material_properties.GetInvSigT()),
    q(material_properties.GetQ()),
    q_per_ster(material_properties.GetQPerSter()),
    is_material_fissile(material_properties.GetFissileIDMap()),
    nu_sigf(material_properties.GetNuSigF()),
    sigs(material_properties.GetSigS()),
    sigs_per_ster(material_properties.GetSigSPerSter()),
    fiss_transfer(material_properties.GetChiNuSigF()),
    fiss_transfer_per_ster(material_properties.GetChiNuSigFPerSter())
{}

template <int dim>
FundamentalData<dim>::FundamentalData (dealii::ParameterHandler &prm,
    dealii::Triangulation<dim> &tria)
    :
    pcout(std::cout,
        (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)),
    //material(prm),
    mesh(prm),
    mat_vec(std::shared_ptr<MatrixVector> (new MatrixVector())),
    //xsec(std::shared_ptr<XSections> (new XSections(material))),
    fe_data(prm),
    dof_handler(tria) {
  MaterialProperties material_properties{prm};
  xsec = std::make_shared<XSections>(material_properties);
  bbuilders::BuildAQ<dim>(prm, aq);
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
FEData<dim>::~FEData() {}

template struct FundamentalData<1>;
template struct FundamentalData<2>;
template struct FundamentalData<3>;

template struct FEData<1>;
template struct FEData<2>;
template struct FEData<3>;

#include "computing_data.h"
#include "bart_builders.h"

#include <deal.II/fe/update_flags.h>

XSections::XSections (std::unique_ptr<MaterialProperties> material)
    :
    sigt(material->GetSigT()),
    inv_sigt(material->GetInvSigT()),
    q(material->GetQ()),
    q_per_ster(material->GetQPerSter()),
    nu_sigf(material->GetNuSigf()),
    sigs(material->GetSigS()),
    sigs_per_ster(material->GetSigSPerSter()),
    fiss_transfer(material->GetFissTransfer()),
    fiss_transfer_per_ster(material->GetFissTransferPerSter()) {}

template <int dim>
FundamentalData<dim>::FundamentalData (dealii::ParameterHandler &prm,
    dealii::Triangulation<dim> &tria)
    :
    pcout(std::cout,
        (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)),
    aq(bbuilders::BuildAQ(prm)),
    material(bbuilders::BuildMaterial(prm)),
    mesh(bbuilders::BuildMesh(prm)),
    mat_vec(std::shared_ptr<MatrixVector> (new MatrixVector())),
    xsec(std::shared_ptr<XSections> (new XSections(material))),
    fe_data(prm),
    dof_handler(tria) {}

template <int dim>
FundamentalData<dim>::~FundamentalData () {}

template <int dim>
FEData<dim>::FEData (const dealii::ParameterHandler &prm)
    :
    fe(bbuilders::BuildFESpaces(prm)) {
  const dealii::UpdateFlags update_flags =
      update_values | update_gradients |
      update_quadrature_points |
      update_JxW_values;

  for (auto &f : fe) {
    // f.first is the equation name, f.second is the correspoding FiniteElement
    // object
    q_rule[f.first] = std::shared_ptr<dealii::QGauss<dim>>(
        new dealii::QGauss<dim>(p_order[f.first]+1));
    qf_rule[f.first] = std::shared_ptr<dealii::QGauss<dim-1>>(
        new dealii::QGauss<dim-1>(p_order[f.first]+1));

    fv[f.first] = std::shared_ptr<dealii::FEValues<dim>> (
        new dealii::FEValues<dim> (*f.second, *q_rule[f.first], update_flags));

    fvf[f.first] = std::shared_ptr<dealii::FEFaceValues<dim>> (
        new dealii::FEFaceValues<dim>(
            *f.second, *qf_rule[f.first], update_flags));

    fvf_nei[f.first] = std::shared_ptr<dealii::FEFaceValues<dim>> (
        new dealii::FEFaceValues<dim>(
            *f.second, *qf_rule[f.first], update_flags));

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
          new dealii::FEValues<dim> (*f.second, *q_rule["corr"], update_flags));

      fvf["corr"] = std::shared_ptr<dealii::FEFaceValues<dim>> (
          new dealii::FEFaceValues<dim>(
              *f.second, *qf_rule["corr"], update_flags));

      fvf_nei["corr"] = std::shared_ptr<dealii::FEFaceValues<dim>> (
          new dealii::FEFaceValues<dim>(
              *f.second, *qf_rule["corr"], update_flags));

      n_q["corr"] = q_rule["corr"]->size();
      n_qf["corr"] = qf_rule["corr"]->size();
    }
  }
}

template struct FundamentalData<1>;
template struct FundamentalData<2>;
template struct FundamentalData<3>;

template struct FEData<1>;
template struct FEData<2>;
template struct FEData<3>;

template struct IterationData<1>;
template struct IterationData<2>;
template struct IterationData<3>;

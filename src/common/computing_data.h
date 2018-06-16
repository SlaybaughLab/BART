#ifndef BART_SRC_COMMON_COMPUTING_DATA_H_
#define BART_SRC_COMMON_COMPUTING_DATA_H_

#include "../material/material_properties.h"
#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"

#include <deal.II/base/parameter_handler.h>

template <int dim>
struct FundamentalData {
  /*!
   Constructor. When constructing, data will be built using builder functions.
   */
  FundamentalData (dealii::ParameterHandler &prm);

  //! Destructor
  ~FundamentalData ();

  dealii::ConditionalOStream pcout;//!< Ostream on one processor.

  std::unique_ptr<AQBase<dim>> aq;//!< Pointer to aq data
  std::unique_ptr<MeshGenerator<dim>> mesh;//!< Pointer to mesh generator
  std::unique_ptr<MaterialProperties> material;//!< Pointer to material properties

  std::shared_ptr<MatrixVector> mat_vec;
  std::shared_ptr<XSections> xsec;//!< Pointer to access all cross sections.

  //! Pointer to finite element related objects of deal::FEValues<dim> etc
  FEData<dim> fe_data;

  //! DoFHandler object
  dealii::DoFHandler<dim> dof_handler;
};

struct MatrixVector {
  std::unordered_map<std::string,
      std::unordered_map<int,
          dealii::PETScWrappers::MPI::SparseMatrix*>> sys_mats;
  std::unordered_map<std::string,
      std::unordered_map<int,
          dealii::PETScWrappers::MPI::Vector*>> sys_flxes;
  std::unordered_map<std::string,
      std::unordered_map<int,
          dealii::PETScWrappers::MPI::Vector*>> sys_rhses;
  std::unordered_map<std::string,
      std::unordered_map<int,
          dealii::PETScWrappers::MPI::Vector*>> sys_fixed_rhses;

  //! Hanging node constraints
  std::unordered_map<std::string,
      std::unordered_map<int, dealii::ConstraintMatrix*>> constraints;

  std::unordered_map<std::string,
      std::map<std::tuple<int,int,int>, dealii::Vector<double>>> moments;
};

struct XSections {
  XSections (std::unique_ptr<MaterialProperties> material);

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> sigt;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> inv_sigt;

  //! \f$Q\f$ values of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> q;

  //! \f$Q/(4\pi)\f$ values of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> q_per_ster;

  //! \f$\nu\sigma_\mathrm{f}\f$ of all groups for all fissile materials.
  /*!
   \todo Change data type to std::unordered_map.
   */
  const std::unordered_map<int, std::vector<double>> nu_sigf;

  //! Scattering matrices for all materials (i.e. \f$\sigma_\mathrm{s,g'\to g}\f$).
  /*!
   \todo Change data type to std::vector<FullMatrix<double> >
   */
  const std::unordered_map<int, dealii::FullMatrix<double>> sigs;

  //! \f$\sigma_\mathrm{s,g'\to g}/(4\pi)\f$
  const std::unordered_map<int, dealii::FullMatrix<double>> sigs_per_ster;

  //! \f$\chi\nu\sigma_\mathrm{f}\f$ of all incident and outgoing groups for fissile materials.
  const std::unordered_map<int, dealii::FullMatrix<double>> fiss_transfer;

  //! \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for fissile materials.
  const std::unordered_map<int, dealii::FullMatrix<double>> fiss_transfer_per_ster;
};

template <int dim>
struct FEData {
  FEData (const dealii::ParameterHandler &prm);
  ~FEData ();

  //!< Finite element spaces.
  std::unordered_map<std::string, dealii::FiniteElement<dim, dim>*> fe;

  // "c" in the following quantities means "correction" for NDA use
  //!< Pointer of quadrature rules in cell.
  std::unordered_map<std::string, std::shared_ptr<dealii::QGauss<dim>>> q_rule;

  //!< Pointer of quadrature rule on cell face.
  std::unordered_map<std::string, std::shared_ptr<dealii::QGauss<dim-1>>> qf_rule;

  //! Pointer of FEValues object.
  /*!
   In short FEValues and FEFaceValues represent finite element evaluated in
   quadrature points of a cell or on a face. For details, please refer to <a href
   ="https://www.dealii.org/8.5.0/doxygen/deal.II/classFEValuesBase.html"
   style="color:blue"><b>FEValues page</b></a>.
   */
  std::unordered_map<std::string, std::shared_ptr<dealii::FEValues<dim>>> fv;

  //! Pointer of FEFaceValues object.
  std::unordered_map<std::string, std::shared_ptr<dealii::FEFaceValues<dim>>> fvf;

  //! Pointer of FEFaceValues object used in DFEM interface bilinear form assembly.
  std::unordered_map<std::string, std::shared_ptr<dealii::FEFaceValues<dim>>> fvf_nei;
};

#endif //BART_SRC_COMMON_COMPUTING_DATA_H_

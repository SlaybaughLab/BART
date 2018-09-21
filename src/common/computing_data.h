#ifndef BART_SRC_COMMON_COMPUTING_DATA_H_
#define BART_SRC_COMMON_COMPUTING_DATA_H_

#include "../material/materials.h"
#include "../material/material_properties_I.h"
#include "../aqdata/aq_base.h"
#include "../mesh/mesh_generator.h"

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/conditional_ostream.h>

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
  XSections (Materials &material);
  XSections (MaterialPropertiesI &material_properties);

  //! \f$\sigma_\mathrm{t}\f$ of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> sigt;

  //! \f$1/\sigma_\mathrm{t}\f$ of all groups for all materials.
  std::unordered_map<int, std::vector<double>> inv_sigt;

  //! \f$Q\f$ values of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> q;

  //! \f$Q/(4\pi)\f$ values of all groups for all materials.
  const std::unordered_map<int, std::vector<double>> q_per_ster;

  const std::unordered_map<int, bool> is_material_fissile;

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

  //! Finite element spaces.
  std::unordered_map<std::string, dealii::FiniteElement<dim, dim>*> fe;

  //! Finite element polynomial orders.
  std::unordered_map<std::string, int> p_order;

  //! Finite element methods
  std::unordered_map<std::string, std::string> discretization;

  // "c" in the following quantities means "correction" for NDA use
  //! Pointer of quadrature rules in cell.
  std::unordered_map<std::string, std::shared_ptr<dealii::QGauss<dim>>> q_rule;

  //! Pointer of quadrature rule on cell face.
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

  std::unordered_map<std::string, int> dofs_per_cell;
  std::unordered_map<std::string, int> n_q;
  std::unordered_map<std::string, int> n_qf;

  std::unordered_map<std::string, std::vector<int>> local_dof_indices;
  std::unordered_map<std::string, std::vector<int>> neigh_dof_indices;
};

template <int dim>
struct FundamentalData {
  /*!
   Constructor. When constructing, data will be built using builder functions.
   */
  FundamentalData (dealii::ParameterHandler &prm,
      dealii::Triangulation<dim> &tria);

  //! Destructor
  ~FundamentalData ();

  dealii::ConditionalOStream pcout;//!< Ostream on one processor.

  std::unique_ptr<AQBase<dim>> aq;//!< Pointer to aq data
  MeshGenerator<dim> mesh;//!< Pointer to mesh generator
  Materials material;//!< Pointer to material properties

  std::shared_ptr<MatrixVector> mat_vec;
  std::shared_ptr<XSections> xsec;//!< Pointer to access all cross sections.

  //! Pointer to finite element related objects of deal::FEValues<dim> etc
  FEData<dim> fe_data;

  std::vector<typename dealii::DoFHandler<dim>::active_cell_iterator> local_cells;

  //! DoFHandler object
  dealii::DoFHandler<dim> dof_handler;
};

#endif //BART_SRC_COMMON_COMPUTING_DATA_H_

#ifndef BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_HPP_
#define BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_HPP_

#include <map>
#include <string>
#include <unordered_map>

#include <deal.II/base/parameter_handler.h>

#include "problem/parameters_i.hpp"
#include "problem/parameter_types.hpp"

namespace bart::problem {

/*!
 * \brief Problem parameters derived using a dealii ParameterHandler object.
 * 
 * The ParameterHandler can be given directly or a file that can be parsed into
 * a ParameterHandler object can be passed.
 */

class ParametersDealiiHandler : public ParametersI {
 public:
  using K_EffectiveUpdaterName = eigenvalue::k_eigenvalue::K_EffectiveUpdaterName;
  struct KeyWords {
    /*!
     * \brief Data struct to contain keywords for input files
     * These are the strings that the ParameterHandler will look for in provided
     * input files.
     */

    // Equation parameters
    const std::string kNEnergyGroups_{ "number of groups" }; // change to energy groups
    const std::string kTransportModel_{ "transport model"} ;

    // Domain parameters
    const std::string kHaveReflectiveBC_{ "have reflective boundary" };
    const std::string kNCells_{ "number of cells for x, y, z directions" };
    const std::string kReflectiveBoundary_{ "reflective boundary names" };
    const std::string kSpatialDimension_{ "problem dimension" };
    const std::string kSpatialMax_{ "x, y, z max values of boundary locations" };
    const std::string kUniformRefinements_{ "uniform refinements" };

    // FEM Parameters
    const std::string kDiscretization_{ "ho spatial discretization" };
    const std::string kFEPolynomialDegree_{ "finite element polynomial degree" };

    // Material parameters
    const std::string kMaterialSubsection_{ "material ID map" };
    const std::string kMaterialMapFilename_{ "material id file name" };
    const std::string kMaterialFilenames_{ "material id file name map" };
    const std::string kNumberOfMaterials_{ "number of materials" };

    // Acceleration parameters
    const std::string kUseTwoGridAcceleration_{ "use two-grid acceleration" };
    const std::string kDoNDA_{ "do nda" };

    // Solver parameters
    const std::string kEigenSolver_{ "eigen solver name" };
    const std::string kK_EffectiveUpdaterType_{ "k_effective updater type" };
    const std::string kInGroupSolver_{ "in group solver name" };
    const std::string kLinearSolver_{ "ho linear solver name" };

    // Quadrature
    const std::string kAngularQuad_{ "angular quadrature name" };
    const std::string kAngularQuadOrder_{ "angular quadrature order" };

    // Instrumentation and output
    const std::string kOutputFilenameBase_{ "output file name base" };
    const std::string kOutputAggregatedSourceData_{ "output aggregated source data"};
    const std::string kOutputScalarFluxAsVTU_{ "output scalar flux as vtu"};
    const std::string kOutputFissionSourceAsVTU_{ "output fission source as vtu"};
    const std::string kOutputScatteringSourceAsVTU_{ "output scattering source as vtu"};
    const std::string kOutputInnerIterationsToFile_{ "output inner iterations to file"};

    // Fourier analysis
    const std::string kDoDFTOfError_{ "do dft of error" };
  };
  
  ParametersDealiiHandler() = default;
  /*! Constructor that parses a given filename in the appropriate format to be
   * read by a ParameterHandler object.
   */
  explicit ParametersDealiiHandler(const std::string& filename);

  /*! \brief Parses a ParameterHandler object for problem parameters */
  auto Parse(dealii::ParameterHandler&) -> void;
  
  /*! \brief Set up a ParameterHandler object for problem parameters.
   * This setup includes default values.
   */
  auto SetUp(dealii::ParameterHandler&) const -> void;

  // Equation parameters
  auto NEnergyGroups() const -> int override { return n_groups_; }
  auto TransportModel() const -> EquationType override { return transport_model_; }

  // Domain parameters
  auto HaveReflectiveBC() const -> bool override { return have_reflective_bc_; }
  auto NCells() const -> std::vector<int> override { return n_cells_; }
  auto ReflectiveBoundary() const -> std::map<Boundary, bool> override { return reflective_boundary_; }
  auto SpatialDimension() const -> int override { return spatial_dimension_; }
  auto SpatialMax() const -> std::vector<double> override { return spatial_max; }
  auto UniformRefinements() const -> int override { return uniform_refinements_; }

  // FEM Parameters
  auto Discretization() const -> DiscretizationType override { return discretization_; }
  auto FEPolynomialDegree() const -> int override { return fe_polynomial_degree_; }

  // Material parameters
  auto MaterialMapFilename() const -> std::string override { return material_map_filename_; }
  auto MaterialFilenames() const -> std::unordered_map<int, std::string> override { return material_filenames_; }
  auto NumberOfMaterials() const -> int override { return n_materials_; }

  // Acceleration parameters
  auto UseTwoGridAcceleration() const -> bool override { return use_two_grid_acceleration_; };
  auto DoNDA() const -> bool override { return do_nda_; }

  // Solver parameters
  auto EigenSolver() const -> EigenSolverType override { return eigen_solver_; }
  auto K_EffectiveUpdaterType() const -> K_EffectiveUpdaterName override { return k_effective_updater_type_; };
  auto InGroupSolver() const -> InGroupSolverType override { return in_group_solver_; }
  auto LinearSolver() const -> LinearSolverType override { return linear_solver_; }

  // Quadrature parameters
  auto AngularQuad() const -> AngularQuadType override { return angular_quad_; }
  auto AngularQuadOrder() const -> int override { return angular_quad_order_; }

  // Instrumentation and output
  auto DoDiscreteFourierTransformOfError() const noexcept -> bool override { return do_dft_of_error_; }
  auto OutputFilenameBase() const noexcept -> std::string override { return output_filename_base_; }
  auto OutputAggregatedSourceData() const -> bool override { return output_aggregated_source_data_; }
  auto OutputScalarFluxAsVTU() const -> bool override { return output_scalar_flux_as_vtu_; }
  auto OutputFissionSourceAsVTU() const -> bool override { return output_fission_source_as_vtu_; }
  auto OutputScatteringSourceAsVTU() const -> bool override { return output_scattering_source_as_vtu_; }
  auto OutputInnerIterationsToFile() const -> bool override { return output_inner_iterations_to_file_; }

  auto GetKeyWords() const noexcept -> KeyWords { return key_words_; }
  
 private:
  // Equation parameters
  int                                  n_groups_{ 0 };
  EquationType                         transport_model_{ EquationType::kNone };
  // Domain parameters
  bool                                 have_reflective_bc_{ false };
  std::vector<int>                     n_cells_{};
  std::map<Boundary, bool>             reflective_boundary_{};
  int                                  spatial_dimension_{ 0 };
  std::vector<double>                  spatial_max{};
  int                                  uniform_refinements_{ 0 };
  // FEM Parameters
  DiscretizationType                   discretization_{ DiscretizationType::kNone };
  int                                  fe_polynomial_degree_{ 0 };
  // Material parameters
  int                                  n_materials_{ 0 };
  std::string                          material_map_filename_{};
  std::unordered_map<int, std::string> material_filenames_{};
  // Acceleration parameters
  bool                                 use_two_grid_acceleration_{ false };
  bool                                 do_nda_{ false };
  // Solver parameters
  EigenSolverType                      eigen_solver_{ EigenSolverType::kNone };
  K_EffectiveUpdaterName               k_effective_updater_type_{ K_EffectiveUpdaterName::kCalculatorViaFissionSource };
  InGroupSolverType                    in_group_solver_{ InGroupSolverType::kNone };
  LinearSolverType                     linear_solver_{ LinearSolverType::kNone };
  //Quadrature
  AngularQuadType                      angular_quad_{ AngularQuadType::kNone };
  int                                  angular_quad_order_{ 0 };

  // Instrumentation and output
  std::string                          output_filename_base_{};
  bool                                 do_dft_of_error_{ false };
  bool                                 output_aggregated_source_data_{ false };
  bool                                 output_scalar_flux_as_vtu_{ false };
  bool                                 output_fission_source_as_vtu_{ false };
  bool                                 output_scattering_source_as_vtu_{ false };
  bool                                 output_inner_iterations_to_file_{ false };

  // Key-words struct                  
  KeyWords                             key_words_{};

  // Options mapping
  // The KeyWords struct provides the "key words" that dealii identifies when looking for options, these are
  // the options and what they map to.

  const std::unordered_map<std::string, Boundary> kBoundaryMap_ {
    {"xmin", Boundary::kXMin}, {"xmax", Boundary::kXMax},
    {"ymin", Boundary::kYMin}, {"ymax", Boundary::kYMax},
    {"zmin", Boundary::kZMin}, {"zmax", Boundary::kZMax},
    }; /*!< Maps boundaries to strings used in parsed input files. */

  const std::unordered_map<std::string, DiscretizationType> kDiscretizationTypeMap_ {
    {"none", DiscretizationType::kNone},
    {"cfem", DiscretizationType::kContinuousFEM},
    {"dfem", DiscretizationType::kDiscontinuousFEM}, }; //!< Maps discretization type to strings used in parsed input


  const std::unordered_map<std::string, EquationType> kEquationTypeMap_ {
      {"none", EquationType::kNone},
      {"diffusion", EquationType::kDiffusion},
      {"saaf",      EquationType::kSelfAdjointAngularFlux} }; /*!< Maps equation type to strings used in parsed input files. */

  const std::unordered_map<std::string, EigenSolverType> kEigenSolverTypeMap_ {
    {"pi",   EigenSolverType::kPowerIteration},
    {"none", EigenSolverType::kNone},
        }; /*!< Maps eigen solver type to strings used in parsed input files. */

  const std::unordered_map<std::string, K_EffectiveUpdaterName> kK_EffectiveUpdaterNameMap_ {
      {"fission source", K_EffectiveUpdaterName::kCalculatorViaFissionSource},
      {"rayleigh quotient", K_EffectiveUpdaterName::kCalculatorViaRayleighQuotient}
  };

  const std::unordered_map<std::string, InGroupSolverType>
  kInGroupSolverTypeMap_ {
    {"si",   InGroupSolverType::kSourceIteration},
    {"none", InGroupSolverType::kNone},
        }; /*!< Maps in-group solver type to strings used in parsed input
            * files. */
  
  const std::unordered_map<std::string, LinearSolverType> kLinearSolverTypeMap_ {
    {"gmres",    LinearSolverType::kGMRES},
    {"none",     LinearSolverType::kNone},
        };  /*!< Maps linear solver type to strings used in parsed input
             * files. */

  const std::unordered_map<std::string, AngularQuadType> kAngularQuadTypeMap_ {
    {"level_symmetric_gaussian", AngularQuadType::kLevelSymmetricGaussian},
    {"gauss_legendre",           AngularQuadType::kGaussLegendre},
    {"none",                     AngularQuadType::kNone},
  }; /*!< Maps angular quadrature type to strings used in parsed input files. */

  // Setup functions
  /*! \brief Set up basic problem parameters */
  auto SetUpBasicParameters(dealii::ParameterHandler&) const -> void;
  /*! \brief Set up mesh parameters */
  auto SetUpMeshParameters(dealii::ParameterHandler&) const -> void;
  /*! \brief Set up material parameters */
  auto SetUpMaterialParameters(dealii::ParameterHandler&) const -> void;
  /*! \brief Set up acceleration parameters */
  auto SetUpAccelerationParameters(dealii::ParameterHandler&) const -> void;
  /*! \brief Set up solver parameters */
  auto SetUpSolverParameters(dealii::ParameterHandler&) const -> void;
  /*! \brief Set up angular quadrature parameters */
  auto SetUpAngularQuadratureParameters(dealii::ParameterHandler&) const -> void;
  
  
  /*! \brief Parses a ParameterHandler entry of type dealii::Patterns::List with
   * doubles into a vector.
   */
  auto ParseDealiiList(std::string to_parse) const -> std::vector<double>;
  auto ParseDealiiIntList(std::string to_parse) const -> std::vector<int>;

  /*! \brief Parses the Material Filename mapping which is of the following
   * deal.II pattern: Map(Integer, Anything)
   */
  auto ParseMap(std::string to_parse) const -> std::unordered_map<int, std::string>;

  /*! \brief Parses a ParameterHandler entry of type
   * dealii::Patterns::MultipleSelection returning a map of one type to another
   */
  template<typename Key>
  auto ParseDealiiMultiple(std::string to_parse,
                           std::unordered_map<std::string, Key> enum_map) const -> std::map<Key, bool>;

  /*! \brief Returns a string formed by combining the key strings in a mapping.
   * 
   * Entries will be separated by `|`. Used to generate valid option strings for
   * ParameterHandler entries. Optional parameter allows a list of options to
   * ignore.
   */
  template<typename T>
  auto GetOptionString(std::unordered_map<std::string, T> enum_map, std::vector<T> to_ignore) const -> std::string;
  template<typename T>
  auto GetOptionString(std::unordered_map<std::string, T> enum_map, T to_ignore) const -> std::string;
  template<typename T>
  auto GetOptionString(std::unordered_map<std::string, T> enum_map) const -> std::string;
};

} // namespace bart::problem

#endif // BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_HPP_

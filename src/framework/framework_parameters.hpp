#ifndef BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_
#define BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_

#include "eigenvalue/k_eigenvalue/factory.hpp"
#include "problem/parameter_types.hpp"
#include "data/cross_sections/material_cross_sections.hpp"
#include "data/cross_sections/one_group_cross_sections_i.hpp"
#include "utility/named_type.h"
#include "quadrature/calculators/angular_flux_integrator_i.hpp"
#include "quadrature/quadrature_set_i.hpp"
#include "system/moments/spherical_harmonic_i.h"
#include "system/solution/solution_types.h"

#include <memory>
#include <optional>
#include <set>

namespace bart::framework {

struct FrameworkParameters {
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorI;
  using AngularFluxStorage = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using AngularQuadratureOrder = quadrature::Order;
  using DomainSize = utility::NamedType<std::vector<double>, struct DomainSizeStruct>;
  using K_EffectiveUpdaterName = eigenvalue::k_eigenvalue::K_EffectiveUpdaterName;
  using NumberOfCells = utility::NamedType<std::vector<int>, struct NumberOfCellsStruct>;
  using PolynomialDegree = utility::NamedType<int, struct PolynomialDegreeStruct>;
  using SpatialDimension = utility::NamedType<int, struct SpatialDimensionStruct>;
  using Moments = system::moments::SphericalHarmonicI;

  std::string name{""};
  std::string output_filename_base{""};

  // Basic properties
  int                         neutron_energy_groups{1};
  problem::EquationType       equation_type{problem::EquationType::kDiffusion};
  std::set<problem::Boundary> reflective_boundaries{};

  // Material properties
  std::string material_mapping{ "" };

  // Solver structure
  std::optional<problem::EigenSolverType> eigen_solver_type{std::nullopt};
  K_EffectiveUpdaterName                  k_effective_updater{ K_EffectiveUpdaterName::kCalculatorViaFissionSource };
  problem::InGroupSolverType              group_solver_type{problem::InGroupSolverType::kSourceIteration};

  // Angular quadrature parameters
  problem::AngularQuadType angular_quadrature_type{ problem::AngularQuadType::kNone };
  std::optional<AngularQuadratureOrder> angular_quadrature_order{ std::nullopt };

  // Solver domain
  SpatialDimension spatial_dimension{ 1 };
  DomainSize domain_size{ {10.0} };
  NumberOfCells number_of_cells { {10} };
  int uniform_refinements { 0 };


  problem::DiscretizationType discretization_type{problem::DiscretizationType::kContinuousFEM};
  problem::CellFiniteElementType cell_finite_element_type{problem::CellFiniteElementType::kGaussian};
  PolynomialDegree polynomial_degree{ 2 };



  // Optional shared framework parts
  template <typename Part> using OptionalSharedPart = std::optional<std::shared_ptr<Part>>;
  OptionalSharedPart<data::cross_sections::CrossSectionsI> cross_sections_ {std::nullopt };

  // Acceleration methods
  bool use_nda_{ false };
  bool use_two_grid_{ false };
  // Indicates "level" of the framework, with 0 being the top level
  int framework_level_{ 0 };
  // Higher order data to support NDA
  struct NDA_Data {
    std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr_{ nullptr };
    std::shared_ptr<Moments> higher_order_moments_ptr_{ nullptr };
    AngularFluxStorage higher_order_angular_flux_{};
  };
  struct TwoGridData {
    std::shared_ptr<data::cross_sections::OneGroupCrossSectionsI> one_group_cross_sections_ptr_{ nullptr };
  };
  NDA_Data nda_data_{};
  TwoGridData two_grid_data_{};

  // Instrumentation options
  bool output_aggregated_source_data{ false };
  bool output_scalar_flux_as_vtu{ false };
  bool output_fission_source_as_vtu{ false };
  bool output_scattering_source_as_vtu{ false };
  bool output_inner_iterations_to_file{ false };
};

} // namespace bart::framework

#endif //BART_SRC_FRAMEWORK_FRAMEWORK_PARAMETERS_HPP_

#ifndef BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_HPP_

#include "calculator/two_grid/spectral_shape_i.hpp"
#include "solver/eigenvalue/spectral_radius_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::calculator::two_grid {

class SpectralShape : public SpectralShapeI, public utility::HasDependencies {
 public:
  using EigenvalueSolver = solver::eigenvalue::SpectralRadiusI;
  using SpectralShapeI::DealiiMatrix; // dealii::FullMatrix<double>;
  SpectralShape(std::unique_ptr<EigenvalueSolver>);
  auto CalculateSpectralShape(const DealiiMatrix &sigma_t, const DealiiMatrix &sigma_s) -> std::vector<double> override;

  auto eigenvalue_solver_ptr() -> EigenvalueSolver*{ return eigenvalue_solver_ptr_.get(); }
 protected:
  std::unique_ptr<EigenvalueSolver> eigenvalue_solver_ptr_{ nullptr };
};

} // namespace bart::calculator::two_grid

#endif //BART_SRC_CALCULATOR_TWO_GRID_SPECTRAL_SHAPE_HPP_

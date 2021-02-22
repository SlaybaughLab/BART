#ifndef BART_SRC_SOLVER_EIGENVALUE_KRYLOV_SCHUR_HPP_
#define BART_SRC_SOLVER_EIGENVALUE_KRYLOV_SCHUR_HPP_

#include "solver/eigenvalue/spectral_radius_i.hpp"

namespace bart::solver::eigenvalue {

class KrylovSchurEigenvalueSolver : public SpectralRadiusI {
 public:
  auto SpectralRadius(const MatrixBase *base) -> std::pair<double, std::vector<double>> override;
};

} // namespace bart::solver::eigenvalue

#endif //BART_SRC_SOLVER_EIGENVALUE_KRYLOV_SCHUR_HPP_

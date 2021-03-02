#ifndef BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_HPP_
#define BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_HPP_

#include "formulation/stamper_i.hpp"
#include "formulation/updater/fixed_updater_i.h"

namespace bart::formulation::updater {

template <int dim>
class FixedUpdater : public FixedUpdaterI {
 public:
  using Matrix = dealii::FullMatrix<double>;
  using Vector = dealii::Vector<double>;
  using CellPtr = domain::CellPtr<dim>;
  using FaceIndex = domain::FaceIndex;
  using MatrixFunction = std::function<void(Matrix&, const CellPtr&)>;
  using VectorFunction = std::function<void(Vector&, const CellPtr&)>;
  using MatrixBoundaryFunction = std::function<void(Matrix&, const FaceIndex, const CellPtr&)>;
  using VectorBoundaryFunction = std::function<void(Vector&, const FaceIndex, const CellPtr&)>;
  using Stamper = formulation::StamperI<dim>;

  FixedUpdater(std::shared_ptr<Stamper> stamper_ptr) : stamper_ptr_(stamper_ptr) {}
  virtual ~FixedUpdater() = default;
  auto UpdateFixedTerms(system::System& to_update, system::EnergyGroup group, quadrature::QuadraturePointIndex index)
  -> void override {
    SetUpFixedFunctions(to_update, group, index);
    const int energy_group{ group.get() };
    auto fixed_matrix_ptr = to_update.left_hand_side_ptr_->GetFixedTermPtr({energy_group, 0});
    auto fixed_vector_ptr = to_update.right_hand_side_ptr_->GetFixedTermPtr({energy_group, 0});
    *fixed_matrix_ptr = 0;
    *fixed_vector_ptr = 0;

    for (auto& matrix_function : fixed_matrix_functions_)
      stamper_ptr_->StampMatrix(*fixed_matrix_ptr, matrix_function);
    for (auto& vector_function : fixed_vector_functions_)
      stamper_ptr_->StampVector(*fixed_vector_ptr, vector_function);
    for (auto& matrix_boundary_function : fixed_matrix_boundary_functions_)
      stamper_ptr_->StampBoundaryMatrix(*fixed_matrix_ptr, matrix_boundary_function);
    // LCOV_EXCL_START
    // Excluded from coverage because no current formulations utilize this.
    for (auto& vector_boundary_function : fixed_vector_boundary_functions_)
      stamper_ptr_->StampBoundaryVector(*fixed_vector_ptr, vector_boundary_function);
    // LCOV_EXCL_STOP
  }

 protected:
  virtual auto SetUpFixedFunctions(system::System&, system::EnergyGroup, quadrature::QuadraturePointIndex) -> void = 0;
  std::vector<MatrixFunction> fixed_matrix_functions_{};
  std::vector<VectorFunction> fixed_vector_functions_{};
  std::vector<MatrixBoundaryFunction> fixed_matrix_boundary_functions_{};
  std::vector<VectorBoundaryFunction> fixed_vector_boundary_functions_{};

  std::shared_ptr<Stamper> stamper_ptr_{ nullptr };

};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_HPP_

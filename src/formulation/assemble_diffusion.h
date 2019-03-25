//#ifndef BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_
//#define BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_
//
//#include <memory>
//
//#include "formulation/equation/diffusion.h"
//
//#include "data/forward_declarations.h"
//#include "data/system_scalar_fluxes.h"
//#include "domain/definition.h"
//#include "problem/parameter_types.h"
//#include "utility/uncopyable.h"
//
//namespace bart {
//
//namespace formulation {
//
//template <int dim>
//class AssembleDiffusion : private utility::Uncopyable {
// public:
//  using GroupNumber = int;
//  using CellMatrix = dealii::FullMatrix<double>;
//  using CellVector = dealii::Vector<double>;
//
//  enum class TermType {
//    kFixed,
//    kVariable,
//  };
//
//  AssembleDiffusion(
//      std::unique_ptr<equation::Diffusion<dim>> equation,
//      std::unique_ptr<domain::Definition<dim>> domain,
//      std::shared_ptr<data::ScalarSystemMatrixPtrs> system_matrix_ptrs,
//      std::shared_ptr<data::ScalarRightHandSidePtrs> right_hand_side_ptrs,
//      std::unique_ptr<data::ScalarRightHandSidePtrs> fixed_right_hand_side_ptrs,
//      std::map<problem::Boundary, bool> reflective_boundary_map);
//  ~AssembleDiffusion() = default;
//
//  void AssembleFixedBilinearTerms(const GroupNumber group);
//  void AssembleFixedLinearTerms(const GroupNumber group);
//
//  void AssembleFissionTerm(const GroupNumber group,
//                           const double k_effective,
//                           const data::ScalarFluxPtrs &scalar_flux_ptrs,
//                           const TermType term_type);
//  void AssembleScatteringSourceTerm(
//      const GroupNumber group,
//      const data::ScalarFluxPtrs &scalar_flux_ptrs,
//      const TermType term_type);
//
//  void ResetRhs(GroupNumber group) {
//    *(right_hand_side_ptrs_->at(group)) =
//        *(fixed_right_hand_side_ptrs_->at(group));
//  }
//
// private:
//
//  data::RightHandSideVector& GetRhs(TermType term_type, GroupNumber group) {
//    if (term_type == TermType::kFixed)
//      return *(right_hand_side_ptrs_->at(group));
//    else
//      return *(fixed_right_hand_side_ptrs_->at(group));
//  }
//
//  // Unique pointers: equation and solver domain
//  std::unique_ptr<equation::Diffusion<dim>> equation_;
//  std::unique_ptr<domain::Definition<dim>> domain_;
//
//  // Shared pointers -> System data to assemble into
//  std::shared_ptr<data::ScalarSystemMatrixPtrs> system_matrix_ptrs_;
//  std::shared_ptr<data::ScalarRightHandSidePtrs> right_hand_side_ptrs_;
//
//  // Internal storage for fixed values
//  std::unique_ptr<data::ScalarRightHandSidePtrs> fixed_right_hand_side_ptrs_;
//
//  // Other
//  std::map<problem::Boundary, bool> reflective_boundary_map_;
//};
//
//} // namespace formulation
//
//} // namespace bart
//
//#endif // BART_SRC_FORMULATION_ASSEMBLE_SCALAR_H_
#ifndef BART_SRC_FORMULATION_RHS_CONSTANT_HPP_
#define BART_SRC_FORMULATION_RHS_CONSTANT_HPP_

#include <memory>

#include "domain/finite_element/finite_element_i.hpp"
#include "formulation/common/rhs_constant_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::formulation::common {

template <int dim>
class RHSConstant : public RHSConstantI<dim>, public utility::HasDependencies {
 public:
  using typename RHSConstantI<dim>::Vector, typename RHSConstantI<dim>::CellPtr;
  using FiniteElement = domain::finite_element::FiniteElementI<dim>;
  virtual ~RHSConstant() = default;
  RHSConstant(std::shared_ptr<FiniteElement> finite_element_ptr) : finite_element_ptr_(finite_element_ptr) {
    this->AssertPointerNotNull(finite_element_ptr_.get(), "finite element", "RHSConstant constructor");
  }

  auto FillCellConstantTerm(Vector& to_fill, const CellPtr&, const Vector& constant_vector) -> void override;

  auto finite_element_ptr() -> FiniteElement* { return finite_element_ptr_.get(); };
 private:
  std::shared_ptr<FiniteElement> finite_element_ptr_{ nullptr };
};

} // namespace bart::formulation::common

#endif //BART_SRC_FORMULATION_RHS_CONSTANT_HPP_

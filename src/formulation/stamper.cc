#include "formulation/stamper.h"

namespace bart {

namespace formulation {

template<int dim>
Stamper<dim>::Stamper(std::shared_ptr<domain::DefinitionI<dim>> domain_ptr)
    : domain_ptr_(domain_ptr) {
  AssertThrow(domain_ptr_ != nullptr,
      dealii::ExcMessage("Error in constructor of formulation::Stamper, "
                         "provided domain::DefinitionI pointer is null"))
}
template<int dim>
void Stamper<dim>::StampMatrix(
    system::MPISparseMatrix& /*to_stamp*/,
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim> &)> /*stamp_function*/) {

}

template class Stamper<1>;
template class Stamper<2>;
template class Stamper<3>;

} // namespace formulation

} // namespace bart

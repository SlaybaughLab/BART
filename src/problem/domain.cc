#include "domain.h"

namespace bart {

namespace problem {

template <int dim>
Domain<dim>::Domain(std::unique_ptr<domain::MeshI<dim>> &mesh,
                    std::unique_ptr<domain::FiniteElementI<dim>> &finite_element)
    : mesh_(std::move(mesh)),
      finite_element_(std::move(finite_element)) {}

// template class Domain<1>(std::unique_ptr<domain::MeshI<1>> &,
//                          std::unique_ptr<domain::FiniteElementI<1>> &);
// template class Domain<2>(std::unique_ptr<domain::MeshI<2>> &,
//                          std::unique_ptr<domain::FiniteElementI<2>> &);
template class Domain<1>;
template class Domain<2>;

} // namespace bart

} // namespace problem

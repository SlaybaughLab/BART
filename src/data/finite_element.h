#ifndef BART_SRC_DATA_FINITE_ELEMENT_H_
#define BART_SRC_DATA_FINITE_ELEMENT_H_

#include <memory>

#include <deal.II/fe/fe_values.h>

#include "../problem/parameter_types.h"

namespace bart {

namespace data {

template <int dim>
class FiniteElement {
 public:
  FiniteElement(bart::problem::DiscretizationType discretization,
                int polynomial_degree);
  ~FiniteElement() = default;

 private:
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_values_;
};

} // namespace data

} // namespace bart 

#endif // BART_SRC_DATA_FINITE_ELEMENT_H_

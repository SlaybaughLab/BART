#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_

#include <memory>

#include "data/cross_sections.h"
#include "domain/finite_element_i.h"
#include "formulation/types.h"
#include "formulation/equation/transport_i.h"

namespace bart {

namespace formulation {

namespace equation {

template <int dim>
class Transport : public TransportI<dim> {
 public:

  using typename TransportI<dim>::CellPtr;

  Transport(EquationType equation_type, DiscretizationType discretization_type)
      : equation_type_(equation_type),
        discretization_type_(discretization_type) {}
  virtual ~Transport() = default;

  Transport& ProvideCrossSections(
      std::shared_ptr<data::CrossSections> cross_sections) override {
    cross_sections_ = cross_sections;
  }

  Transport& ProvideFiniteElement(
      std::shared_ptr<domain::FiniteElementI<dim>> finite_element) override {
    finite_element_ = finite_element;
    cell_degrees_of_freedom_ = finite_element_->dofs_per_cell();
    cell_quadrature_points_ = finite_element_->n_cell_quad_pts();
    return *this;
  }

  EquationType equation_type() const override { return equation_type_; };
  DiscretizationType discretization_type() const override {
    return discretization_type_; };

  Transport& SetCell(const CellPtr &to_set) {
    if (finite_element_->values()->get_cell() != to_set)
      finite_element_->values()->reinit(to_set);
    return *this;
  }

 protected:
  EquationType equation_type_;
  DiscretizationType discretization_type_;
  std::shared_ptr<domain::FiniteElementI<dim>> finite_element_;
  std::shared_ptr<data::CrossSections> cross_sections_;

  int cell_degrees_of_freedom_ = 0;
  int cell_quadrature_points_ = 0;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
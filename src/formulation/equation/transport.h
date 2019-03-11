#ifndef BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
#define BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_

#include <memory>

#include "data/system_scalar_fluxes.h"
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
  using typename TransportI<dim>::FaceNumber;

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
    face_quadrature_points_ = finite_element_->n_face_quad_pts();
    return *this;
  }

  Transport& ProvideScalarFluxes(
      std::shared_ptr<data::ScalarFluxPtrs> scalar_fluxes) override {
    scalar_fluxes_ = scalar_fluxes;
    return *this;
  }

  EquationType equation_type() const override { return equation_type_; };
  DiscretizationType discretization_type() const override {
    return discretization_type_; };

  void SetCell(const CellPtr &to_set) const override {
    if (finite_element_->values()->get_cell() != to_set)
      finite_element_->values()->reinit(to_set);
  }

  void SetFace(const CellPtr &to_set, FaceNumber face) const override {
    if (finite_element_->face_values()->get_cell() != to_set &&
        finite_element_->face_values()->get_face_index() != face)
      finite_element_->face_values()->reinit(to_set, face);
  }

 protected:
  EquationType equation_type_;
  DiscretizationType discretization_type_;
  mutable std::shared_ptr<domain::FiniteElementI<dim>> finite_element_;
  std::shared_ptr<data::CrossSections> cross_sections_;
  std::shared_ptr<data::ScalarFluxPtrs> scalar_fluxes_;

  int cell_degrees_of_freedom_ = 0;
  int cell_quadrature_points_ = 0;
  int face_quadrature_points_ = 0;
};

} // namespace equation

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_EQUATION_TRANSPORT_H_
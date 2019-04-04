#include "domain/finite_element.h"

namespace bart {

namespace domain {

template <int dim>
bool FiniteElement<dim>::SetCell(const CellPtr &to_set) {
  bool already_set = false;

  if (values_reinit_called_) {
    already_set = (values_->get_cell()->id() == to_set->id());
  }

  if (!already_set) {
    values_->reinit(to_set);
    values_reinit_called_ = true;
  }

  return !already_set;
}

template <int dim>
bool FiniteElement<dim>::SetFace(const CellPtr &to_set, const FaceNumber face) {
  bool already_set = false;
  bool cell_already_set = false;

  if (face_values_reinit_called_) {
    cell_already_set = (face_values_->get_cell()->id() == to_set->id());
    bool face_already_set =
        (static_cast<int>(face_values_->get_face_index()) == face);

    already_set = (cell_already_set && face_already_set);
  }

  if (!already_set) {
    face_values_->reinit(to_set, face);
    face_values_reinit_called_ = true;
  }

  return !cell_already_set;
}
template<int dim>
std::vector<double> FiniteElement<dim>::ValueAtQuadrature(
    const data::MomentVector moment) const {

  std::vector<double> return_vector(finite_element_->dofs_per_cell, 0);

  values_->get_function_values(moment, return_vector);

  return return_vector;
}

template class FiniteElement<1>;
template class FiniteElement<2>;
template class FiniteElement<3>;

} // namespace domain

} // namespace bart
#include "domain/finite_element.h"

namespace bart {

namespace domain {

template <int dim>
bool FiniteElement<dim>::SetCell(const CellPtr &to_set) {
  bool already_set = false;

  if (reinit_has_been_called_) {
  already_set = (values()->get_cell()->id() == to_set->id());
  }

  if (!already_set) {
  values()->reinit(to_set);
  reinit_has_been_called_ = true;
  }

  return !already_set;
}

template <int dim>
bool FiniteElement<dim>::SetFace(const CellPtr &to_set, const FaceNumber face) {
  bool already_set = false;
  bool cell_already_set = false;

  if (reinit_has_been_called_) {
  cell_already_set = (face_values()->get_cell()->id() == to_set->id());
  bool face_already_set =
      (static_cast<int>(face_values()->get_face_index()) == face);

  already_set = (cell_already_set && face_already_set);
  }

  if (!already_set) {
  face_values()->reinit(to_set, face);
  reinit_has_been_called_ = true;
  }

  return !cell_already_set;
}

template class FiniteElement<1>;
template class FiniteElement<2>;
template class FiniteElement<3>;

} // namespace domain

} // namespace bart
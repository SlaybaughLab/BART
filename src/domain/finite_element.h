#ifndef BART_DOMAIN_FINITE_ELEMENT_H_
#define BART_DOMAIN_FINITE_ELEMENT_H_

#include "domain/finite_element_i.h"

namespace bart {

namespace domain {

template<int dim>
class FiniteElement : public FiniteElementI<dim> {
 public:
  using typename FiniteElementI<dim>::CellPtr;
  using typename FiniteElementI<dim>::FaceNumber;

  virtual ~FiniteElement() = default;

  bool SetCell(const CellPtr &to_set) override {
    bool already_set = (values()->get_cell() == to_set);

    if (!already_set)
      values()->reinit(to_set);

    return already_set;
  }

  bool SetFace(const FaceNumber face) override {

    bool already_set =
        (static_cast<int>(face_values()->get_face_index()) == face);

    if (!already_set) {
      auto cell = values()->get_cell();
      face_values()->reinit(cell, face);
    }

    return already_set;
  }

  bool SetFace(const CellPtr &to_set, const FaceNumber face) override {
    bool cell_already_set = (values()->get_cell() == to_set);
    bool face_already_set =
        (static_cast<int>(face_values()->get_face_index()) == face);
    bool already_set = (!cell_already_set && !face_already_set);

    if (!already_set) {
      face_values()->reinit(to_set, face);
    }

    return already_set;
  }

 protected:
  using FiniteElementI<dim>::values;
  using FiniteElementI<dim>::face_values;

};


} // namespace domain

} // namespace bart

#endif // BART_DOMAIN_FINITE_ELEMENT_H_
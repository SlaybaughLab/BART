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
    return false;
  }
//    if (finite_element_->face_values()->get_cell() != to_set &&
//        static_cast<int>(finite_element_->face_values()->get_face_index())
//            != face)
//      finite_element_->face_values()->reinit(to_set, face);
//  }

 protected:
  using FiniteElementI<dim>::values;

};


} // namespace domain

} // namespace bart

#endif // BART_DOMAIN_FINITE_ELEMENT_H_
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

  // Basic Finite Element data
  dealii::FiniteElement<dim, dim> *finite_element() override {
    return finite_element_.get(); };

  int dofs_per_cell() const override { return finite_element_->dofs_per_cell; }

  int n_cell_quad_pts() const override { return cell_quadrature_->size(); }

  int n_face_quad_pts() const override { return face_quadrature_->size(); }

  dealii::FEValues<dim> *values() override { return values_.get(); };

  dealii::FEFaceValues<dim> *face_values() override {
    return face_values_.get(); };

  dealii::FEFaceValues<dim> *neighbor_face_values() override {
    return neighbor_face_values_.get(); };

  dealii::QGauss<dim> *cell_quadrature() { return cell_quadrature_.get(); };

  dealii::QGauss<dim - 1> *face_quadrature() { return face_quadrature_.get(); };


  bool SetCell(const CellPtr &to_set) override {
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

  bool SetFace(const CellPtr &to_set, const FaceNumber face) override {
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

  double ShapeValue(const int cell_degree_of_freedom,
                    const int cell_quadrature_point) const override {
    return 0;
  }

 protected:
  std::shared_ptr<dealii::FiniteElement<dim, dim>> finite_element_;
  std::shared_ptr<dealii::FEValues<dim>> values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> face_values_;
  std::shared_ptr<dealii::FEFaceValues<dim>> neighbor_face_values_;
  std::shared_ptr<dealii::QGauss<dim>> cell_quadrature_;
  std::shared_ptr<dealii::QGauss<dim - 1>> face_quadrature_;

  bool reinit_has_been_called_ = false;
  using FiniteElementI<dim>::values;
  using FiniteElementI<dim>::face_values;

};


} // namespace domain

} // namespace bart

#endif // BART_DOMAIN_FINITE_ELEMENT_H_
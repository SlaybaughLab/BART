#include "formulation/equation/diffusion.h"

namespace bart {

namespace formulation {

namespace equation {

template<int dim>
void Diffusion<dim>::Precalculate() {

}

template<int dim>
void Diffusion<dim>::FillCellBilinearTerm(Matrix &to_fill,
                                          const CellPtr &cell_ptr,
                                          const GroupNumber group) const {
  MaterialID material_id = cell_ptr->material_id();

  double sigma_t = cross_sections_->sigma_t.at(material_id)[group];
  double sigma_s = cross_sections_->sigma_s.at(material_id)(group, group);
  double sigma_r = sigma_t - sigma_s;

  double diff_coef = cross_sections_->diffusion_coef.at(material_id)[group];

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i,j) += 0;
      }
    }
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;

} // namespace equation

} // namespace formulation

} // namespace bart
//#include "formulation/equation/diffusion.h"
//#include "diffusion.h"
//
//#include <vector>
//
//namespace bart {
//
//namespace formulation {
//
//namespace equation {
//
//template<int dim>
//void Diffusion<dim>::Precalculate(const CellPtr &cell_ptr) {
//  SetCell(cell_ptr);
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    FullMatrix shape_squared(cell_degrees_of_freedom_,
//                             cell_degrees_of_freedom_);
//    FullMatrix gradient_squared(cell_degrees_of_freedom_,
//                                cell_degrees_of_freedom_);
//
//    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
//      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
//        shape_squared(i,j) = finite_element_->values()->shape_value(i,q) *
//            finite_element_->values()->shape_value(j,q);
//        gradient_squared(i,j) = finite_element_->values()->shape_grad(i,q) *
//            finite_element_->values()->shape_grad(j,q);
//      }
//    }
//    shape_squared_[q] = shape_squared;
//    gradient_squared_[q] = gradient_squared;
//  }
//}
//
//template<int dim>
//void Diffusion<dim>::FillCellStreamingTerm(FullMatrix &to_fill,
//                                           const CellPtr &cell_ptr,
//                                           const GroupNumber group) const {
//  SetCell(cell_ptr);
//  MaterialID material_id = cell_ptr->material_id();
//
//  double diff_coef = cross_sections_->diffusion_coef.at(material_id)[group];
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    double gradient = finite_element_->values()->JxW(q);
//    to_fill.add(diff_coef*gradient, gradient_squared_[q]);
//  }
//}
//
//template <int dim>
//void Diffusion<dim>::FillCellCollisionTerm(FullMatrix &to_fill,
//                                           const CellPtr &cell_ptr,
//                                           const GroupNumber group) const {
//  SetCell(cell_ptr);
//  MaterialID material_id = cell_ptr->material_id();
//
//  double sigma_t = cross_sections_->sigma_t.at(material_id)[group];
//  double sigma_s = cross_sections_->sigma_s.at(material_id)(group, group);
//  double sigma_r = sigma_t - sigma_s;
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    double gradient = finite_element_->values()->JxW(q);
//    to_fill.add(sigma_r*gradient, shape_squared_[q]);
//  }
//}
//
//template<int dim>
//void Diffusion<dim>::FillBoundaryTerm(
//    FullMatrix &to_fill,
//    const CellPtr &cell_ptr,
//    const GroupNumber,
//    const FaceNumber face_number,
//    const BoundaryType boundary_type) const {
//
//  if (boundary_type == BoundaryType::kVacuum) {
//    SetFace(cell_ptr, face_number);
//
//    for (int q = 0; q < face_quadrature_points_; ++q) {
//      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
//        for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
//          to_fill(i,j) += (
//              finite_element_->face_values()->shape_value(i, q) *
//              finite_element_->face_values()->shape_value(j,q)
//              ) * 0.5 * finite_element_->face_values()->JxW(q);
//        }
//      }
//    }
//  }
//}
//
//template<int dim>
//void Diffusion<dim>::FillCellFixedSource(Vector &rhs_to_fill,
//                                         const CellPtr &cell_ptr,
//                                         const GroupNumber group) const {
//  SetCell(cell_ptr);
//  MaterialID material_id = cell_ptr->material_id();
//
//  // Fixed source at each cell quadrature point
//  std::vector<double> cell_fixed_source(cell_quadrature_points_);
//
//  cell_fixed_source = std::vector<double>(
//      cell_quadrature_points_,
//      cross_sections_->q.at(material_id)[group]);
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    cell_fixed_source[q] *= finite_element_->values()->JxW(q);
//    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
//      rhs_to_fill(i) += finite_element_->values()->shape_value(i, q) *
//          cell_fixed_source[q];
//    }
//  }
//}
//template<int dim>
//void Diffusion<dim>::FillCellScatteringSource(
//    Vector &rhs_to_fill,
//    const CellPtr &cell_ptr,
//    const GroupNumber group,
//    const data::MomentsMap& out_group_moments) const {
//
//  SetCell(cell_ptr);
//  MaterialID material_id = cell_ptr->material_id();
//
//  std::vector<double> cell_scatter_source(cell_quadrature_points_);
//
//  for (auto const &moment_pair : out_group_moments) {
//    auto &[moment_index, moment] = moment_pair;
//    int group_in = moment_index[0];
//
//    if ((group_in != group) && (moment_index[1] == 0)
//        && (moment_index[2] == 0)) {
//
//      std::vector<double> group_cell_scalar_flux(cell_quadrature_points_);
//      double
//          sigma_s = cross_sections_->sigma_s.at(material_id)(group_in, group);
//
//      if (sigma_s >= 1e-14) {
//        finite_element_->values()->get_function_values(moment,
//                                                       group_cell_scalar_flux);
//
//        for (int q = 0; q < cell_quadrature_points_; ++q) {
//          cell_scatter_source[q] += sigma_s * group_cell_scalar_flux[q];
//
//        }
//      }
//    }
//  }
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    cell_scatter_source[q] *= finite_element_->values()->JxW(q);
//    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
//      rhs_to_fill(i) += finite_element_->values()->shape_value(i, q) *
//          cell_scatter_source[q];
//    }
//  }
//}
//
//template <int dim>
//void Diffusion<dim>::FillCellFissionSource(
//    Vector &rhs_to_fill,
//    const CellPtr &cell_ptr,
//    const GroupNumber group,
//    const double k_effective,
//    const data::MomentsMap &group_moments) const {
//
//  SetCell(cell_ptr);
//  MaterialID material_id = cell_ptr->material_id();
//
//  std::vector<double> cell_fission_source(cell_quadrature_points_);
//  std::vector<double> group_cell_scalar_flux(cell_quadrature_points_);
//
//  finite_element_->values()->get_function_values(group_moments[{group, 0, 0}],
//                                                 group_cell_scalar_flux);
//
//  double scaled_fission_transfer =
//      cross_sections_->fiss_transfer.at(material_id)(group, group)/k_effective;
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    cell_fission_source[q] += scaled_fission_transfer * group_cell_scalar_flux[q];
//  }
//
//  for (int q = 0; q < cell_quadrature_points_; ++q) {
//    cell_fission_source[q] *= finite_element_->values()->JxW(q);
//    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
//      rhs_to_fill(i) += finite_element_->values()->shape_value(i, q) *
//          cell_fission_source[q];
//    }
//  }
//}
//
//template class Diffusion<1>;
//template class Diffusion<2>;
//template class Diffusion<3>;
//
//} // namespace equation
//
//} // namespace formulation
//
//} // namespace bart
#include "data/cross_sections/collapsed_one_group_cross_sections.hpp"

#include <numeric>

namespace bart::data::cross_sections {

namespace  {
template <typename MappedType> using MaterialIDMappedTo = std::unordered_map<int, MappedType>;
using FullMatrix = dealii::FullMatrix<double>;

auto collapse_vector = [](MaterialIDMappedTo<std::vector<double>> to_collapse) {
  MaterialIDMappedTo<std::vector<double>> return_map;
  for (auto& [id, vector] : to_collapse)
    return_map[id] = std::vector<double>(1, std::accumulate(vector.cbegin(), vector.cend(), 0.0));
  return return_map;
};

auto collapse_matrix = [](MaterialIDMappedTo<FullMatrix> to_collapse) {
  MaterialIDMappedTo<FullMatrix> return_map;
  for (auto& [id, matrix] : to_collapse) {
    return_map[id] = FullMatrix(1, 1);
    for (const auto val : matrix)
      return_map[id](0, 0) += val;
  }
  return return_map;
};

} // namespace

CollapsedOneGroupCrossSections::CollapsedOneGroupCrossSections(const CrossSectionsI &to_collapse) {
  diffusion_coef_ = collapse_vector(to_collapse.diffusion_coef());
  sigma_t_ = collapse_vector(to_collapse.sigma_t());
  inverse_sigma_t_ = collapse_vector(to_collapse.inverse_sigma_t());
  sigma_s_ = collapse_matrix(to_collapse.sigma_s());
  sigma_s_per_ster_ = collapse_matrix(to_collapse.sigma_s_per_ster());
  q_ = collapse_vector(to_collapse.q());
  q_per_ster_ = collapse_vector(to_collapse.q_per_ster());
  is_material_fissile_ = to_collapse.is_material_fissile();
  nu_sigma_f_ = collapse_vector(to_collapse.nu_sigma_f());
  fiss_transfer_ = collapse_matrix(to_collapse.fiss_transfer());
  fiss_transfer_per_ster_ = collapse_matrix(to_collapse.fiss_transfer_per_ster());
}

} // namespace bart::data::cross_sections

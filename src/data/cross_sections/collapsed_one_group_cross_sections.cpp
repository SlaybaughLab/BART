#include "data/cross_sections/collapsed_one_group_cross_sections.hpp"

#include <numeric>

namespace bart::data::cross_sections {

namespace  {
template <typename MappedType> using MaterialIDMappedTo = std::unordered_map<int, MappedType>;
using FullMatrix = dealii::FullMatrix<double>;

auto ScaleVector = [](MaterialIDMappedTo<std::vector<double>> to_scale, MaterialIDMappedTo<std::vector<double>> scaling) {
  for (auto& [id, vector] : to_scale) {
    for (std::vector<double>::size_type group = 0; group < vector.size(); ++group) {
      vector.at(group) *= scaling.at(id).at(group);
    }
  }
  return to_scale;
};

auto collapse_vector = [](MaterialIDMappedTo<std::vector<double>> to_collapse) {
  MaterialIDMappedTo<std::vector<double>> return_map;
  for (auto& [id, vector] : to_collapse)
    return_map[id] = std::vector<double>(1, std::accumulate(vector.cbegin(), vector.cend(), 0.0));
  return return_map;
};

auto ScaleMatrix = [](MaterialIDMappedTo<FullMatrix> to_scale, MaterialIDMappedTo<std::vector<double>> scaling) {
  for (auto& [id, matrix] : to_scale) {
    for (FullMatrix::size_type group = 0; group < matrix.m(); ++group) {
      for (FullMatrix::size_type group_in = 0; group_in < matrix.n(); ++group_in) {
        matrix(group, group_in) *= scaling.at(id).at(group_in);
      }
    }
  }
  return to_scale;
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
  const auto sigma_s{ to_collapse.sigma_s() };
  diffusion_coef_ = collapse_vector(to_collapse.diffusion_coef());
  sigma_t_ = collapse_vector(to_collapse.sigma_t());
  inverse_sigma_t_ = collapse_vector(to_collapse.inverse_sigma_t());
  sigma_s_ = collapse_matrix(sigma_s);
  sigma_s_per_ster_ = collapse_matrix(to_collapse.sigma_s_per_ster());
  q_ = collapse_vector(to_collapse.q());
  q_per_ster_ = collapse_vector(to_collapse.q_per_ster());
  is_material_fissile_ = to_collapse.is_material_fissile();
  nu_sigma_f_ = collapse_vector(to_collapse.nu_sigma_f());
  fiss_transfer_ = collapse_matrix(to_collapse.fiss_transfer());
  fiss_transfer_per_ster_ = collapse_matrix(to_collapse.fiss_transfer_per_ster());

  const int total_groups = sigma_s.begin()->second.m();
  for (const auto& [material_id, sigma_t_vector] : sigma_t_) {
    double sigma_absorption{ sigma_t_vector.at(0) };
    double sigma_removal{ sigma_absorption };
    for (int group = 0; group < total_groups; ++group) {
      sigma_removal -= sigma_s.at(material_id)(group, group);
      for (int group_in = 0; group_in < total_groups; ++group_in) {
        sigma_absorption -= sigma_s.at(material_id)(group, group_in);
      }
    }
    sigma_absorption_[material_id] = sigma_absorption;
    sigma_removal_[material_id] = sigma_removal;
  }
}
CollapsedOneGroupCrossSections::CollapsedOneGroupCrossSections(
    const CrossSectionsI &to_collapse,
    const MaterialIDMappedTo<std::vector<double>>& scaling_factor_by_group) {
  const auto sigma_s{ to_collapse.sigma_s() };
  diffusion_coef_ = collapse_vector(ScaleVector(to_collapse.diffusion_coef(), scaling_factor_by_group));
  sigma_t_ = collapse_vector(ScaleVector(to_collapse.sigma_t(), scaling_factor_by_group));
  auto inverse_scaling_factor_by_group{ scaling_factor_by_group };
  for (auto& pair : inverse_scaling_factor_by_group) {
    for (auto& value : pair.second)
      value = 1.0/value;
  }
  //diffusion_coef_ = collapse_vector(ScaleVector(to_collapse.diffusion_coef(), inverse_scaling_factor_by_group));
  inverse_sigma_t_ = collapse_vector(ScaleVector(to_collapse.inverse_sigma_t(), inverse_scaling_factor_by_group));
  sigma_s_ = collapse_matrix(ScaleMatrix(sigma_s, scaling_factor_by_group));
  sigma_s_per_ster_ = collapse_matrix(ScaleMatrix(to_collapse.sigma_s_per_ster(), scaling_factor_by_group));
  q_ = collapse_vector(ScaleVector(to_collapse.q(), scaling_factor_by_group));
  q_per_ster_ = collapse_vector(ScaleVector(to_collapse.q_per_ster(), scaling_factor_by_group));
  is_material_fissile_ = to_collapse.is_material_fissile();
  nu_sigma_f_ = collapse_vector(ScaleVector(to_collapse.nu_sigma_f(), scaling_factor_by_group));

  MaterialIDMappedTo<FullMatrix> transposed_fiss_transfer;
  MaterialIDMappedTo<FullMatrix> transposed_fiss_transfer_per_ster;

  const auto fission_transfer = to_collapse.fiss_transfer();
  const auto fission_transfer_per_ster = to_collapse.fiss_transfer_per_ster();

  for (const auto& [id, matrix] : fission_transfer) {
    transposed_fiss_transfer[id] = FullMatrix();
    transposed_fiss_transfer.at(id).copy_transposed(matrix);
    transposed_fiss_transfer_per_ster[id] = FullMatrix();
    transposed_fiss_transfer_per_ster.at(id).copy_transposed(fission_transfer_per_ster.at(id));
  }

  fiss_transfer_ = collapse_matrix(ScaleMatrix(transposed_fiss_transfer, scaling_factor_by_group));
  fiss_transfer_per_ster_ = collapse_matrix(ScaleMatrix(transposed_fiss_transfer_per_ster, scaling_factor_by_group));

  const int total_groups = sigma_s.begin()->second.m();
  for (const auto& [material_id, sigma_t_vector] : sigma_t_) {
    double sigma_absorption{ sigma_t_vector.at(0) };
    double sigma_removal{ sigma_absorption };
    for (int group = 0; group < total_groups; ++group) {
      const double within_group_scattering_{ sigma_s.at(material_id)(group, group) * scaling_factor_by_group.at(material_id).at(group) };
      sigma_removal -= within_group_scattering_;
      for (int group_in = 0; group_in < total_groups; ++group_in) {
        sigma_absorption -= sigma_s.at(material_id)(group, group_in) * scaling_factor_by_group.at(material_id).at(group_in);
      }
    }
    sigma_absorption_[material_id] = sigma_absorption;
    sigma_removal_[material_id] = sigma_removal;
  }
}

} // namespace bart::data::cross_sections

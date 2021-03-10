#include "data/cross_sections/collapsed_one_group_cross_sections.hpp"

#include <numeric>

#include "data/cross_sections/tests/cross_sections_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;
using ::testing::NiceMock, ::testing::Return, ::testing::DoDefault, ::testing::ContainerEq;

class DataCrossSectionsCollapsedOneGroup : public ::testing::Test {
 public:
  using CrossSectionsMock = NiceMock<data::cross_sections::CrossSectionsMock>;
  using FullMatrix = dealii::FullMatrix<double>;
  template <typename MappedType> using MaterialIDMappedTo = std::unordered_map<int, MappedType>;

  std::shared_ptr<CrossSectionsMock> cross_sections_mock_ptr_{ std::make_shared<CrossSectionsMock>() };

  // Test parameters
  static constexpr int total_groups{ 3 };
  static constexpr int total_materials{ 2 };

  MaterialIDMappedTo<std::vector<double>> collapsed_diffusion_coef_, scaled_collapsed_diffusion_coef_;
  MaterialIDMappedTo<std::vector<double>> collapsed_sigma_t_, scaled_collapsed_sigma_t_;
  MaterialIDMappedTo<std::vector<double>> collapsed_inverse_sigma_t_, scaled_collapsed_inverse_sigma_t_;
  MaterialIDMappedTo<FullMatrix> collapsed_sigma_s_, scaled_collapsed_sigma_s_;
  MaterialIDMappedTo<FullMatrix> collapsed_sigma_s_per_ster_, scaled_collapsed_sigma_s_per_ster_;
  MaterialIDMappedTo<std::vector<double>> collapsed_q_, scaled_collapsed_q_;
  MaterialIDMappedTo<std::vector<double>> collapsed_q_per_ster_, scaled_collapsed_q_per_ster_;
  MaterialIDMappedTo<bool> is_material_fissile_;
  MaterialIDMappedTo<std::vector<double>> collapsed_nu_sigma_f_, scaled_collapsed_nu_sigma_f_;
  MaterialIDMappedTo<FullMatrix> collapsed_fiss_transfer_, scaled_collapsed_fiss_transfer_;
  MaterialIDMappedTo<FullMatrix> collapsed_fiss_transfer_per_ster_, scaled_collapsed_fiss_transfer_per_ster_;

  MaterialIDMappedTo<double> sigma_r_, scaled_sigma_r_;
  MaterialIDMappedTo<double> sigma_a_, scaled_sigma_a_;


  MaterialIDMappedTo<std::vector<double>> scaling_factor_by_group_, inverse_scaling_factor_by_group_;


  auto SetUp() -> void override;
};

auto DataCrossSectionsCollapsedOneGroup::SetUp() -> void {
  auto ScaleAndCollapseVector = [total_groups = this->total_groups](MaterialIDMappedTo<std::vector<double>> to_collapse,
                                                                    const MaterialIDMappedTo<std::vector<double>> scaling) {
    MaterialIDMappedTo<std::vector<double>> return_map;
    for (auto& [id, vector] : to_collapse) {
      double sum{ 0.0 };
      for (int group = 0; group < total_groups; ++group)
        sum += vector.at(group) * scaling.at(id).at(group);
      return_map[id] = std::vector<double>(1, sum);
    }
    return return_map;
  };
  auto CollapseVector = [=, total_groups = this->total_groups](MaterialIDMappedTo<std::vector<double>> to_collapse) {
    MaterialIDMappedTo<std::vector<double>> scaling_factor{{0, {1, 1, 1}}, {1, {1, 1, 1}}};
    return ScaleAndCollapseVector(to_collapse, scaling_factor);
  };
  auto ScaleAndCollapseMatrix = [G = this->total_groups](MaterialIDMappedTo<FullMatrix> to_collapse,
                                                         const MaterialIDMappedTo<std::vector<double>> scaling) {
    MaterialIDMappedTo<FullMatrix> return_map;
    for (auto& [id, matrix] : to_collapse) {
      FullMatrix material_matrix(1, 1);

      for (int g = 0; g < G; ++g) {
        for (int g_in = 0; g_in < G; ++g_in) {
          material_matrix(0, 0) += matrix(g, g_in) * scaling.at(id).at(g_in);
        }
      }
      return_map[id] = material_matrix;
    }
    return return_map;
  };
  auto CollapseMatrix = [=, G = this->total_groups](MaterialIDMappedTo<FullMatrix> to_collapse) {
    MaterialIDMappedTo<std::vector<double>> scaling_factor{{0, {1, 1, 1}}, {1, {1, 1, 1}}};
    return ScaleAndCollapseMatrix(to_collapse, scaling_factor);
  };
  auto RandomMaterialVectorMap = [total_materials = this->total_materials, total_groups = this->total_groups]() {
    MaterialIDMappedTo<std::vector<double>> return_map;
    for (int material_id = 0; material_id < total_materials; ++material_id) {
      return_map[material_id] = test_helpers::RandomVector(total_groups, 0, 100);
    }
    return return_map;
  };
  auto RandomMaterialMatrixMap = [total_materials = this->total_materials, total_groups = this->total_groups]() {
    MaterialIDMappedTo<FullMatrix> return_map;
    for (int material_id = 0; material_id < total_materials; ++material_id) {
      return_map[material_id] = test_helpers::RandomMatrix(total_groups, total_groups, 0, 100);
    }
    return return_map;
  };

  MaterialIDMappedTo<std::vector<double>> sigma_t{{0, {0.86, 0.83, 0.05}}, {1, {0.82, 0.1, 0.79}}};
  MaterialIDMappedTo<std::vector<std::vector<double>>> sigma_s_values{
      {0, {{0.51, 0.96, 0.10}, {0.26, 0.67, 0.20}, {0.16, 0.98, 0.43}}},
      {1, {{0.96, 0.44, 0.94}, {0.58, 0.88, 0.40}, {0.27, 0.28, 0.42}}}
  };
  MaterialIDMappedTo<FullMatrix> sigma_s;
  for (int mat_id = 0; mat_id < total_materials; ++mat_id) {
    is_material_fissile_[mat_id] = test_helpers::RandomInt(0, 1) > 0.5;
    sigma_s[mat_id] = FullMatrix(total_groups, total_groups);
    for (int m = 0; m < total_groups; ++m) {
      for (int n = 0; n < total_groups; ++n) {
        sigma_s.at(mat_id).set(m, n, sigma_s_values.at(mat_id).at(m).at(n));
      }
    }
  }
  sigma_r_ = {{0, 0.13}, {1, -0.55}};
  sigma_a_ = {{0, -1.13}, {1, -2.33}};

  scaling_factor_by_group_ = {{0, {0.22, 0.11, 0.67}}, {1, {0.07, 0.79, 0.14}}};
  inverse_scaling_factor_by_group_ = {{0, {1.0/0.22, 1.0/0.11, 1.0/0.67}}, {1, {1.0/0.07, 1.0/0.79, 1.0/0.14}}};

  scaled_collapsed_sigma_t_ = {{0, {0.314}}, {1, {0.247}}};
  scaled_sigma_r_ = {{0, -0.16}, {1, -0.5742}};
  scaled_sigma_a_ = {{0, -0.4666}, {1, -1.1094}};

  auto diffusion_coef = RandomMaterialVectorMap();
  collapsed_diffusion_coef_ = CollapseVector(diffusion_coef);
  scaled_collapsed_diffusion_coef_ = ScaleAndCollapseVector(diffusion_coef, scaling_factor_by_group_);

  collapsed_sigma_t_ = CollapseVector(sigma_t);

  auto inverse_sigma_t = RandomMaterialVectorMap();
  collapsed_inverse_sigma_t_ = CollapseVector(inverse_sigma_t);
  scaled_collapsed_inverse_sigma_t_ = ScaleAndCollapseVector(inverse_sigma_t, inverse_scaling_factor_by_group_);

  collapsed_sigma_s_ = CollapseMatrix(sigma_s);
  scaled_collapsed_sigma_s_ = ScaleAndCollapseMatrix(sigma_s, scaling_factor_by_group_);

  const auto sigma_s_per_ster = RandomMaterialMatrixMap();
  collapsed_sigma_s_per_ster_ = CollapseMatrix(sigma_s_per_ster);
  scaled_collapsed_sigma_s_per_ster_ = ScaleAndCollapseMatrix(sigma_s_per_ster, scaling_factor_by_group_);

  auto q = RandomMaterialVectorMap();
  collapsed_q_ = CollapseVector(q);
  scaled_collapsed_q_ = ScaleAndCollapseVector(q, scaling_factor_by_group_);

  auto q_per_ster = RandomMaterialVectorMap();
  collapsed_q_per_ster_ = CollapseVector(q_per_ster);
  scaled_collapsed_q_per_ster_ = ScaleAndCollapseVector(q_per_ster, scaling_factor_by_group_);

  auto nu_sigma_f = RandomMaterialVectorMap();
  collapsed_nu_sigma_f_ = CollapseVector(nu_sigma_f);
  scaled_collapsed_nu_sigma_f_ = ScaleAndCollapseVector(nu_sigma_f, scaling_factor_by_group_);

  auto fission_transfer = RandomMaterialMatrixMap();
  collapsed_fiss_transfer_ = CollapseMatrix(fission_transfer);
  scaled_collapsed_fiss_transfer_ = ScaleAndCollapseMatrix(fission_transfer, scaling_factor_by_group_);

  auto fission_transfer_per_ster = RandomMaterialMatrixMap();
  collapsed_fiss_transfer_per_ster_ = CollapseMatrix(fission_transfer_per_ster);
  scaled_collapsed_fiss_transfer_per_ster_ = ScaleAndCollapseMatrix(fission_transfer_per_ster, scaling_factor_by_group_);

  ON_CALL(*cross_sections_mock_ptr_, diffusion_coef()).WillByDefault(Return(diffusion_coef));
  ON_CALL(*cross_sections_mock_ptr_, sigma_t()).WillByDefault(Return(sigma_t));
  ON_CALL(*cross_sections_mock_ptr_, inverse_sigma_t()).WillByDefault(Return(inverse_sigma_t));
  ON_CALL(*cross_sections_mock_ptr_, sigma_s()).WillByDefault(Return(sigma_s));
  ON_CALL(*cross_sections_mock_ptr_, sigma_s_per_ster()).WillByDefault(Return(sigma_s_per_ster));
  ON_CALL(*cross_sections_mock_ptr_, q()).WillByDefault(Return(q));
  ON_CALL(*cross_sections_mock_ptr_, q_per_ster()).WillByDefault(Return(q_per_ster));
  ON_CALL(*cross_sections_mock_ptr_, is_material_fissile()).WillByDefault(Return(is_material_fissile_));
  ON_CALL(*cross_sections_mock_ptr_, nu_sigma_f()).WillByDefault(Return(nu_sigma_f));
  ON_CALL(*cross_sections_mock_ptr_, fiss_transfer()).WillByDefault(Return(fission_transfer));
  ON_CALL(*cross_sections_mock_ptr_, fiss_transfer_per_ster()).WillByDefault(Return(fission_transfer_per_ster));
}

TEST_F(DataCrossSectionsCollapsedOneGroup, Constructor) {
  EXPECT_CALL(*cross_sections_mock_ptr_, diffusion_coef()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_t()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, inverse_sigma_t()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_s()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_s_per_ster()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, q()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, q_per_ster()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, is_material_fissile()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, nu_sigma_f()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, fiss_transfer()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, fiss_transfer_per_ster()).WillOnce(DoDefault());

  data::cross_sections::CollapsedOneGroupCrossSections collapsed_cross_sections(*cross_sections_mock_ptr_);

  EXPECT_THAT(collapsed_cross_sections.diffusion_coef(), ContainerEq(collapsed_diffusion_coef_));
  EXPECT_THAT(collapsed_cross_sections.sigma_t(), ContainerEq(collapsed_sigma_t_));
  EXPECT_THAT(collapsed_cross_sections.inverse_sigma_t(), ContainerEq(collapsed_inverse_sigma_t_));
  EXPECT_THAT(collapsed_cross_sections.sigma_s(), ContainerEq(collapsed_sigma_s_));
  EXPECT_THAT(collapsed_cross_sections.sigma_s_per_ster(), ContainerEq(collapsed_sigma_s_per_ster_));
  EXPECT_THAT(collapsed_cross_sections.q(), ContainerEq(collapsed_q_));
  EXPECT_THAT(collapsed_cross_sections.q_per_ster(), ContainerEq(collapsed_q_per_ster_));
  EXPECT_THAT(collapsed_cross_sections.is_material_fissile(), ContainerEq(is_material_fissile_));
  EXPECT_THAT(collapsed_cross_sections.nu_sigma_f(), ContainerEq(collapsed_nu_sigma_f_));
  EXPECT_THAT(collapsed_cross_sections.fiss_transfer(), ContainerEq(collapsed_fiss_transfer_));
  EXPECT_THAT(collapsed_cross_sections.fiss_transfer_per_ster(), ContainerEq(collapsed_fiss_transfer_per_ster_));
  for (int i = 0; i < total_materials; ++i) {
    EXPECT_NEAR(collapsed_cross_sections.SigmaRemoval().at(i), sigma_r_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaAbsorption().at(i), sigma_a_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaAbsorption(i), sigma_a_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaRemoval(i), sigma_r_.at(i), 1e-10);
  }
}

TEST_F(DataCrossSectionsCollapsedOneGroup, ConstructorWithScaling) {
  EXPECT_CALL(*cross_sections_mock_ptr_, diffusion_coef()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_t()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, inverse_sigma_t()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_s()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, sigma_s_per_ster()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, q()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, q_per_ster()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, is_material_fissile()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, nu_sigma_f()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, fiss_transfer()).WillOnce(DoDefault());
  EXPECT_CALL(*cross_sections_mock_ptr_, fiss_transfer_per_ster()).WillOnce(DoDefault());

  data::cross_sections::CollapsedOneGroupCrossSections collapsed_cross_sections(*cross_sections_mock_ptr_,
                                                                                scaling_factor_by_group_);

  for (int mat_id = 0; mat_id < total_materials; ++mat_id) {
    EXPECT_NEAR(collapsed_cross_sections.diffusion_coef().at(mat_id).at(0), scaled_collapsed_diffusion_coef_.at(mat_id).at(0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.sigma_t().at(mat_id).at(0), scaled_collapsed_sigma_t_.at(mat_id).at(0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.inverse_sigma_t().at(mat_id).at(0), scaled_collapsed_inverse_sigma_t_.at(mat_id).at(0), 1e-10);

    EXPECT_NEAR(collapsed_cross_sections.q().at(mat_id).at(0), scaled_collapsed_q_.at(mat_id).at(0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.q_per_ster().at(mat_id).at(0), scaled_collapsed_q_per_ster_.at(mat_id).at(0), 1e-10);
    EXPECT_EQ(collapsed_cross_sections.is_material_fissile().at(mat_id), is_material_fissile_.at(mat_id));
    EXPECT_NEAR(collapsed_cross_sections.nu_sigma_f().at(mat_id).at(0), scaled_collapsed_nu_sigma_f_.at(mat_id).at(0), 1e-10);

    EXPECT_NEAR(collapsed_cross_sections.sigma_s().at(mat_id)(0,0), scaled_collapsed_sigma_s_.at(mat_id)(0,0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.sigma_s_per_ster().at(mat_id)(0,0), scaled_collapsed_sigma_s_per_ster_.at(mat_id)(0,0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.fiss_transfer().at(mat_id)(0,0), scaled_collapsed_fiss_transfer_.at(mat_id)(0,0), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.fiss_transfer_per_ster().at(mat_id)(0,0), scaled_collapsed_fiss_transfer_per_ster_.at(mat_id)(0,0), 1e-10);
  }
  for (int i = 0; i < total_materials; ++i) {
    EXPECT_NEAR(collapsed_cross_sections.SigmaRemoval().at(i), scaled_sigma_r_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaAbsorption().at(i), scaled_sigma_a_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaAbsorption(i), scaled_sigma_a_.at(i), 1e-10);
    EXPECT_NEAR(collapsed_cross_sections.SigmaRemoval(i), scaled_sigma_r_.at(i), 1e-10);
  }
}

} // namespace

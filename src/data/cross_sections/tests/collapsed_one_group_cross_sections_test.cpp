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
  const int total_groups{ test_helpers::RandomInt(5, 10) };
  const int total_materials{ test_helpers::RandomInt(3, 5) };

  MaterialIDMappedTo<std::vector<double>> collapsed_diffusion_coef_;
  MaterialIDMappedTo<std::vector<double>> collapsed_sigma_t_;
  MaterialIDMappedTo<std::vector<double>> collapsed_inverse_sigma_t_;
  MaterialIDMappedTo<FullMatrix> collapsed_sigma_s_;
  MaterialIDMappedTo<FullMatrix> collapsed_sigma_s_per_ster_;
  MaterialIDMappedTo<std::vector<double>> collapsed_q_;
  MaterialIDMappedTo<std::vector<double>> collapsed_q_per_ster_;
  MaterialIDMappedTo<bool> is_material_fissile_;
  MaterialIDMappedTo<std::vector<double>> collapsed_nu_sigma_f_;
  MaterialIDMappedTo<FullMatrix> collapsed_fiss_transfer_;
  MaterialIDMappedTo<FullMatrix> collapsed_fiss_transfer_per_ster_;

  auto SetUp() -> void override;
};

auto DataCrossSectionsCollapsedOneGroup::SetUp() -> void {
  auto collapse_vector = [](MaterialIDMappedTo<std::vector<double>> to_collapse) {
    MaterialIDMappedTo<std::vector<double>> return_map;
    for (auto& [id, vector] : to_collapse) {
      double sum{ 0.0 };
      for (const double val : vector)
        sum += val;
      return_map[id] = std::vector<double>(1, sum);
    }
    return return_map;
  };

  auto collapse_matrix = [G = this->total_groups](MaterialIDMappedTo<FullMatrix> to_collapse) {
    MaterialIDMappedTo<FullMatrix> return_map;
    for (auto& [id, matrix] : to_collapse) {
      FullMatrix material_matrix(1, 1);

      for (int g = 0; g < G; ++g) {
        for (int g_in = 0; g_in < G; ++g_in) {
          material_matrix(0, 0) += matrix(g, g_in);
        }
      }
      return_map[id] = material_matrix;
    }
    return return_map;
  };

  auto diffusion_coef = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_diffusion_coef_ = collapse_vector(diffusion_coef);
  auto sigma_t = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_sigma_t_ = collapse_vector(sigma_t);
  auto inverse_sigma_t = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_inverse_sigma_t_ = collapse_vector(inverse_sigma_t);
  const auto sigma_s = test_helpers::RandomIntMatrixMap(total_materials, total_groups, total_groups);
  collapsed_sigma_s_ = collapse_matrix(sigma_s);
  const auto sigma_s_per_ster = test_helpers::RandomIntMatrixMap(total_materials, total_groups, total_groups);
  collapsed_sigma_s_per_ster_ = collapse_matrix(sigma_s_per_ster);
  auto q = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_q_ = collapse_vector(q);
  auto q_per_ster = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_q_per_ster_ = collapse_vector(q_per_ster);
  for (int i = 0; i < total_materials; ++i)
    is_material_fissile_[i] = test_helpers::RandomDouble(0, 1) > 0.5;
  auto nu_sigma_f = test_helpers::RandomIntVectorMap(total_materials, total_groups);
  collapsed_nu_sigma_f_ = collapse_vector(nu_sigma_f);
  auto fission_transfer = test_helpers::RandomIntMatrixMap(total_materials, total_groups, total_groups);
  collapsed_fiss_transfer_ = collapse_matrix(fission_transfer);
  auto fission_transfer_per_ster = test_helpers::RandomIntMatrixMap(total_materials, total_groups, total_groups);
  collapsed_fiss_transfer_per_ster_ = collapse_matrix(fission_transfer_per_ster);


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
}

} // namespace

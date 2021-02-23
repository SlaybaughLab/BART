#ifndef BART_SRC_MATERIAL_TESTS_MATERIAL_MOCK_HPP_
#define BART_SRC_MATERIAL_TESTS_MATERIAL_MOCK_HPP_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "data/material/material_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::material {

class MaterialMock : public MaterialI {
 public:
  using int_bool_map = std::unordered_map<int, bool>;
  using int_vector_map = std::unordered_map<int, std::vector<double>>;
  using int_matrix_map = std::unordered_map<int, dealii::FullMatrix<double>>;
  
  MOCK_METHOD(int_bool_map, GetFissileIDMap, (), (const, override));
  MOCK_METHOD(int_vector_map, GetDiffusionCoef, (), (const, override));
  MOCK_METHOD(int_vector_map, GetSigT, (), (const, override));
  MOCK_METHOD(int_vector_map, GetInvSigT, (), (const, override));
  MOCK_METHOD(int_vector_map, GetQ, (), (const, override));
  MOCK_METHOD(int_vector_map, GetQPerSter, (), (const, override));
  MOCK_METHOD(int_vector_map, GetNuSigF, (), (const, override));
  MOCK_METHOD(int_matrix_map, GetSigS, (), (const, override));
  MOCK_METHOD(int_matrix_map, GetSigSPerSter, (), (const, override));
  MOCK_METHOD(int_matrix_map, GetChiNuSigF, (), (const, override));
  MOCK_METHOD(int_matrix_map, GetChiNuSigFPerSter, (), (const, override));
};

} // namespace bart::material

#endif //BART_SRC_MATERIAL_TESTS_MOCK_MATERIAL_PROPERTIES_H_

#ifndef BART_SRC_MATERIAL_TESTS_MATERIAL_MOCK_HPP_
#define BART_SRC_MATERIAL_TESTS_MATERIAL_MOCK_HPP_

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>
#include "../../test_helpers/gmock_wrapper.h"

#include "../material_i.hpp"

namespace btest {

class MockMaterial : public MaterialBase {
 public:

  using int_bool_map = std::unordered_map<int, bool>;
  using int_vector_map = std::unordered_map<int, std::vector<double>>;
  using int_matrix_map = std::unordered_map<int, dealii::FullMatrix<double>>;
  
  MOCK_CONST_METHOD0(GetFissileIDMap, int_bool_map());

  MOCK_CONST_METHOD0(GetDiffusionCoef, int_vector_map());
  
  //! Returns all \f$\sigma_\mathrm{t}\f$ for all groups.
  MOCK_CONST_METHOD0(GetSigT, int_vector_map());

  //! Returns all \f$1/\sigma_\mathrm{t}\f$ for all groups.
  MOCK_CONST_METHOD0(GetInvSigT, int_vector_map());

  //! Returns all fixed source value \f$Q\f$'s for all groups.
  MOCK_CONST_METHOD0(GetQ, int_vector_map());

  //! Returns all \f$Q/(4\pi)\f$'s for all groups.
  MOCK_CONST_METHOD0(GetQPerSter, int_vector_map());

  //! Returns all \f$\nu\sigma_\mathrm{f}\f$'s.
  MOCK_CONST_METHOD0(GetNuSigF, int_vector_map());

  //! Returns all scattering transfer matrices.
  MOCK_CONST_METHOD0(GetSigS, int_matrix_map());

  //! Returns all scattering transfer matrices scaled by \f$4\pi\f$.
  MOCK_CONST_METHOD0(GetSigSPerSter, int_matrix_map());

  //! Returns \f$\chi\nu\sigma_\mathrm{f}\f$ for all fissile materials.
  MOCK_CONST_METHOD0(GetChiNuSigF, int_matrix_map());

  //! Returns \f$\chi\nu\sigma_\mathrm{f}/(4\pi)\f$ for all fissile materials.
  MOCK_CONST_METHOD0(GetChiNuSigFPerSter, int_matrix_map());
  
};

} // namespace btest

#endif //BART_SRC_MATERIAL_TESTS_MOCK_MATERIAL_PROPERTIES_H_

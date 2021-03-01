#ifndef BART_SRC_DATA_CROSS_SECTIONS_TESTS_CROSS_SECTIONS_MOCK_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_TESTS_CROSS_SECTIONS_MOCK_HPP_

#include "data/cross_sections/cross_sections_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::data::cross_sections {

class CrossSectionsMock : public CrossSectionsI {
 public:
  using CrossSectionsI::DealiiMatrix;
  using CrossSectionsI::MaterialIDMappedTo;

  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, diffusion_coef, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, sigma_t, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, inverse_sigma_t, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<DealiiMatrix>, sigma_s, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<DealiiMatrix>, sigma_s_per_ster, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, q, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, q_per_ster, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<bool>, is_material_fissile, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<std::vector<double>>, nu_sigma_f, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<DealiiMatrix>, fiss_transfer, (), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<DealiiMatrix>, fiss_transfer_per_ster, (), (const, override));
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_TESTS_CROSS_SECTIONS_MOCK_HPP_

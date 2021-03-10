#ifndef BART_SRC_DATA_CROSS_SECTIONS_TESTS_ONE_GROUP_CROSS_SECTIONS_MOCK_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_TESTS_ONE_GROUP_CROSS_SECTIONS_MOCK_HPP_

#include "data/cross_sections/one_group_cross_sections_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::data::cross_sections {

class OneGroupCrossSectionsMock : public OneGroupCrossSectionsI {
 public:
  template <typename MappedType>
  using MaterialIDMappedTo = std::unordered_map<int, MappedType>;

  MOCK_METHOD(MaterialIDMappedTo<double>, SigmaRemoval, (), (const, override));
  MOCK_METHOD(double, SigmaRemoval, (int material_id), (const, override));
  MOCK_METHOD(MaterialIDMappedTo<double>, SigmaAbsorption, (), (const, override));
  MOCK_METHOD(double, SigmaAbsorption, (int material_id), (const, override));
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_TESTS_ONE_GROUP_CROSS_SECTIONS_MOCK_HPP_

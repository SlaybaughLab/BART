#ifndef BART_SRC_DATA_CROSS_SECTIONS_ONE_GROUP_CROSS_SECTIONS_I_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_ONE_GROUP_CROSS_SECTIONS_I_HPP_

#include <unordered_map>

namespace bart::data::cross_sections {

class OneGroupCrossSectionsI {
 public:
  template <typename MappedType>
  using MaterialIDMappedTo = std::unordered_map<int, MappedType>;

  ~OneGroupCrossSectionsI() = default;
  virtual auto SigmaRemoval() const -> MaterialIDMappedTo<double> = 0;
  virtual auto SigmaRemoval(int material_id) const -> double = 0;
  virtual auto SigmaAbsorption() const -> MaterialIDMappedTo<double> = 0;
  virtual auto SigmaAbsorption(int material_id) const -> double = 0;
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_ONE_GROUP_CROSS_SECTIONS_I_HPP_

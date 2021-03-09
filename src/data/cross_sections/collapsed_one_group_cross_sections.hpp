#ifndef BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_

#include "data/cross_sections/cross_sections.hpp"
#include "data/cross_sections/one_group_cross_sections_i.hpp"
#include "data/cross_sections/scalable_i.hpp"

namespace bart::data::cross_sections {

class CollapsedOneGroupCrossSections : public CrossSections, public OneGroupCrossSectionsI, public ScalableI {
 public:
  template <typename MappedType>
  using MaterialIDMappedTo = std::unordered_map<int, MappedType>;

  CollapsedOneGroupCrossSections(const CrossSectionsI& to_collapse);
  auto SigmaRemoval() const -> MaterialIDMappedTo<double> override { return sigma_removal_; };
  auto SigmaRemoval(int material_id) const -> double override { return sigma_removal_.at(material_id); };
  auto SigmaAbsorption() const -> MaterialIDMappedTo<double> override { return sigma_absorption_; };
  auto SigmaAbsorption(int material_id) const -> double override { return sigma_absorption_.at(material_id); };

  auto Scale(const MaterialIDMappedTo<std::vector<double>>& scaling_factor_by_group) -> void override;

 protected:
  MaterialIDMappedTo<double> sigma_absorption_, sigma_removal_;
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_

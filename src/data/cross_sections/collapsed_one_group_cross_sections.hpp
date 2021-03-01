#ifndef BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_

#include "data/cross_sections/cross_sections.hpp"

namespace bart::data::cross_sections {

class CollapsedOneGroupCrossSections : public CrossSections {
 public:
  CollapsedOneGroupCrossSections(const CrossSectionsI& to_collapse);
};

} // namespace bart::data::cross_sections

#endif //BART_SRC_DATA_CROSS_SECTIONS_COLLAPSED_ONE_GROUP_CROSS_SECTIONS_HPP_

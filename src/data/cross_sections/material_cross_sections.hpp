#ifndef BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_
#define BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_

#include "cross_sections.hpp"

#include <unordered_map>
#include <vector>

#include <deal.II/lac/full_matrix.h>

#include "data/material/material_i.hpp"

namespace bart::data::cross_sections {

/*! \brief Cross-sections generated from a material.
 *
 * This is the default implementation of cross-sections that parses a MaterialI object for data.
 *
 */
class MaterialCrossSections: public CrossSections {
 public:
  MaterialCrossSections(data::material::MaterialI &materials);
};

} // namespace bart::data::cross_sections

#endif // BART_SRC_DATA_CROSS_SECTIONS_MATERIAL_CROSS_SECTIONS_HPP_

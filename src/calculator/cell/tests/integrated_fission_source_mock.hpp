#ifndef BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_HPP_
#define BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_HPP_

#include "calculator/cell/integrated_fission_source_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::cell {

template<int dim>
class IntegratedFissionSourceMock : public IntegratedFissionSourceI<dim> {
 public:
  MOCK_METHOD(double, CellValue, (domain::CellPtr<dim>, system::moments::SphericalHarmonicI*), (const, override));
};

} // namespace bart::calculator::cell

#endif // BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_HPP_
#ifndef BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_H_
#define BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_H_

#include "calculator/cell/integrated_fission_source_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace calculator {

namespace cell {

template<int dim>
class IntegratedFissionSourceMock : public IntegratedFissionSourceI<dim> {
 public:
  MOCK_CONST_METHOD2_T(CellValue, double(domain::CellPtr<dim>,
      system::moments::SphericalHarmonicI*));

};

} // namespace cell

} // namespace calculator

} // namespace bart

#endif // BART_SRC_CALCULATOR_CELL_TESTS_INTEGRATED_FISSION_SOURCE_MOCK_H_
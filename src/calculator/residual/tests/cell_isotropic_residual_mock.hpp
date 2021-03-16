#ifndef BART_SRC_CALCULATOR_RESIDUAL_TESTS_CELL_ISOTROPIC_RESIDUAL_MOCK_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_TESTS_CELL_ISOTROPIC_RESIDUAL_MOCK_HPP_

#include "calculator/residual/cell_isotropic_residual_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::residual {

template <int dim>
class CellIsotropicResidualMock : public CellIsotropicResidualI<dim> {
 public:
  using typename CellIsotropicResidualI<dim>::CellPtr;
  using typename CellIsotropicResidualI<dim>::FluxMoments;
  MOCK_METHOD(double, CalculateCellResidual, (CellPtr, FluxMoments*, FluxMoments*, int), (override));
};

} // namespace bart::calculator::residual


#endif //BART_SRC_CALCULATOR_RESIDUAL_TESTS_CELL_ISOTROPIC_RESIDUAL_MOCK_HPP_

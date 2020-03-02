#ifndef BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_
#define BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_

#include <functional>
#include <memory>

#include "domain/domain_types.h"
#include "formulation/formulation_types.h"
#include "formulation/tests/stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace formulation {

namespace updater {

namespace test_helpers {

using ::testing::WithArg, ::testing::Invoke, ::testing::_;

template <int dim>
class UpdaterTests {
 public:
  using StamperType = formulation::StamperMock<dim>;

  std::unique_ptr<StamperMock<dim>> MakeStamper();

  static void EvaluateFunction(std::function<void(formulation::FullMatrix&,
                                                  const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateVectorFunction(std::function<void(formulation::Vector&,
                                                        const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateBoundaryFunction(std::function<void(formulation::FullMatrix&,
                                                          const domain::FaceIndex,
                                                          const domain::CellPtr<dim>&)> stamp_function);
  static void EvaluateVectorBoundaryFunction(std::function<void(formulation::Vector&,
                                                                const domain::FaceIndex,
                                                                const domain::CellPtr<dim>&)> stamp_function);
};

template <int dim>
std::unique_ptr<StamperMock<dim>> UpdaterTests<dim>::MakeStamper() {
  auto mock_stamper_ptr = std::make_unique<StamperMock<dim>>();

  ON_CALL(*mock_stamper_ptr, StampMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateFunction)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryMatrix(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateBoundaryFunction)));
  ON_CALL(*mock_stamper_ptr, StampVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateVectorFunction)));
  ON_CALL(*mock_stamper_ptr, StampBoundaryVector(_,_))
      .WillByDefault(WithArg<1>(Invoke(this->EvaluateVectorBoundaryFunction)));

  return std::move(mock_stamper_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  stamp_function(to_stamp, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorFunction(
    std::function<void(formulation::Vector&,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  domain::CellPtr<dim> cell_ptr;
  stamp_function(to_stamp, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateBoundaryFunction(
    std::function<void(formulation::FullMatrix&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::FullMatrix to_stamp;
  domain::CellPtr<dim> cell_ptr;
  domain::FaceIndex index(0);
  stamp_function(to_stamp, index, cell_ptr);
}

template <int dim>
void UpdaterTests<dim>::EvaluateVectorBoundaryFunction(
    std::function<void(formulation::Vector&,
                       const domain::FaceIndex,
                       const domain::CellPtr<dim>&)> stamp_function) {
  formulation::Vector to_stamp;
  domain::CellPtr<dim> cell_ptr;
  domain::FaceIndex index(0);
  stamp_function(to_stamp, index, cell_ptr);
}

} // namespace test_helpers

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_TESTS_UPDATER_TESTS_H_

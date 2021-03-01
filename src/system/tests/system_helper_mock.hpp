#ifndef BART_SRC_SYSTEM_TESTS_SYSTEM_HELPER_MOCK_HPP_
#define BART_SRC_SYSTEM_TESTS_SYSTEM_HELPER_MOCK_HPP_

#include "test_helpers/gmock_wrapper.h"
#include "system/system_helper_i.hpp"

namespace bart::system {

template <int dim>
class SystemHelperMock : public SystemHelperI<dim> {
 public:
  MOCK_METHOD(void, InitializeSystem, (system::System&, const int, const int, const bool, const bool), (const, override));
  MOCK_METHOD(void, SetUpEnergyGroupToAngularSolutionPtrMap, (solution::EnergyGroupToAngularSolutionPtrMap&,
      const int, const int), (const, override));
  MOCK_METHOD(void, SetUpMPIAngularSolution, (system::solution::MPIGroupAngularSolutionI&,
      const domain::DomainI<dim> &, const double), (const, override));
  MOCK_METHOD(void, SetUpSystemTerms, (system::System&, const domain::DomainI<dim>&), (const, override));
  MOCK_METHOD(void, SetUpSystemMoments, (system::System&, const std::size_t), (const, override));
};

} // namespace bart::system

#endif //BART_SRC_SYSTEM_TESTS_SYSTEM_HELPER_MOCK_HPP_

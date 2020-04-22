#ifndef BART_SRC_QUADRATURE_TESTS_QUADRATURE_SET_MOCK_H_
#define BART_SRC_QUADRATURE_TESTS_QUADRATURE_SET_MOCK_H_

#include "quadrature/quadrature_set_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

template <int dim>
class QuadratureSetMock : public QuadratureSetI<dim> {
 public:
  using typename QuadratureSetI<dim>::Iterator;
  using typename QuadratureSetI<dim>::ConstIterator;

  MOCK_METHOD1_T(AddPoint, bool(std::shared_ptr<QuadraturePointI<dim>>));
  MOCK_METHOD2_T(SetReflection, void(std::shared_ptr<QuadraturePointI<dim>>,
      std::shared_ptr<QuadraturePointI<dim>>));
  MOCK_CONST_METHOD1_T(GetReflection, std::shared_ptr<QuadraturePointI<dim>>(
      std::shared_ptr<QuadraturePointI<dim>>));
  MOCK_CONST_METHOD1_T(GetReflectionIndex, std::optional<int> (
      std::shared_ptr<QuadraturePointI<dim>>));
  MOCK_CONST_METHOD1_T(GetQuadraturePoint,
      std::shared_ptr<QuadraturePointI<dim>>(QuadraturePointIndex));
  MOCK_CONST_METHOD1_T(GetQuadraturePointIndex, int(
      std::shared_ptr<QuadraturePointI<dim>>));
  MOCK_CONST_METHOD0_T(quadrature_point_indices, std::set<int>());
  MOCK_METHOD0_T(begin, Iterator());
  MOCK_METHOD0_T(end, Iterator());
  MOCK_CONST_METHOD0_T(cbegin, ConstIterator());
  MOCK_CONST_METHOD0_T(cend, ConstIterator());

  MOCK_CONST_METHOD0_T(size, std::size_t());


};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_TESTS_QUADRATURE_SET_MOCK_H_

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

  MOCK_METHOD(bool, AddPoint, (std::shared_ptr<QuadraturePointI<dim>>), (override));
  MOCK_METHOD(void, SetReflection, (std::shared_ptr<QuadraturePointI<dim>>,
      std::shared_ptr<QuadraturePointI<dim>>), (override));
  MOCK_METHOD(std::shared_ptr<QuadraturePointI<dim>>, GetReflection,
              (std::shared_ptr<QuadraturePointI<dim>>), (override, const));
  MOCK_METHOD(std::optional<int>, GetReflectionIndex,
              (std::shared_ptr<QuadraturePointI<dim>>), (override, const));
  MOCK_METHOD(std::shared_ptr<QuadraturePointI<dim>>, GetQuadraturePoint,
              (QuadraturePointIndex), (override, const));
  MOCK_METHOD(int, GetQuadraturePointIndex,
              (std::shared_ptr<QuadraturePointI<dim>>), (override, const));
  MOCK_METHOD(std::set<int>, quadrature_point_indices, (), (override, const));
  MOCK_METHOD(Iterator, begin, (), (override));
  MOCK_METHOD(Iterator, end, (), (override));
  MOCK_METHOD(ConstIterator, cbegin, (), (override, const));
  MOCK_METHOD(ConstIterator, cend, (), (override, const));
  MOCK_METHOD(std::size_t, size, (), (override, const));


};

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_TESTS_QUADRATURE_SET_MOCK_H_

#include "../aq_base.h"

#include "gtest/gtest.h"

template <int dim>
class AQBaseTest : public AQBase<dim>
{
public:
  AQBaseTest (dealii::ParameterHandler &prm);
  ~AQBaseTest ();
};

template <int dim>
AQBaseTest<dim>::AQBaseTest (dealii::ParameterHandler &prm)
    : AQBase<dim> (prm)
{}

template <int dim>
AQBaseTest<dim>::~AQBaseTest ()
{}

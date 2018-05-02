#include "computing_data.h"
#include "bart_builders.h"

template <int dim>
FundamentalData<dim>::FundamentalData (dealii::ParameterHandler &prm) {
  bbuilders::BuildAQ (prm, aq);
  bbuilders::BuildMaterial (prm, material);
  bbuilders::BuildMesh (prm, mesh);
}

template <int dim>
FundamentalData<dim>::~FundamentalData () {}

template struct FundamentalData<1>;
template struct FundamentalData<2>;
template struct FundamentalData<3>;

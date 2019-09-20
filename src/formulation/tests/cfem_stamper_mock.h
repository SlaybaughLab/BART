#ifndef BART_SRC_FORMULATION_TESTS_CFEM_STAMPER_MOCK_H_
#define BART_SRC_FORMULATION_TESTS_CFEM_STAMPER_MOCK_H_

#include "formulation/cfem_stamper_i.h"

#include "test_helpers/gmock_wrapper.h"


namespace bart {

namespace formulation {

class CFEM_StamperMock : public CFEMStamperI {
 public:
  MOCK_METHOD2(StampStreamingTerm, void(MPISparseMatrix&, const GroupNumber));
  MOCK_METHOD2(StampCollisionTerm, void(MPISparseMatrix&, const GroupNumber));
  MOCK_METHOD1(StampBoundaryTerm, void(MPISparseMatrix&));
  MOCK_METHOD2(StampFixedSource, void(MPIVector&, const GroupNumber));
  MOCK_METHOD5(StampFissionSource,
               void(MPIVector&, const GroupNumber, const double,
                   const system::moments::MomentVector&, const system::moments::MomentsMap&));
  MOCK_METHOD3(StampScatteringSource,
               void(MPIVector&, const GroupNumber,
                   const system::moments::MomentsMap&));
  MOCK_METHOD1(AddReflectiveBoundary, CFEM_StamperMock&(Boundary));
  MOCK_METHOD1(RemoveReflectiveBoundary, CFEM_StamperMock&(Boundary));
  MOCK_CONST_METHOD0(reflective_boundaries, std::unordered_set<Boundary>());
};

} // namespace formulation

} // namespace bart

#endif // BART_SRC_FORMULATION_TESTS_CFEM_STAMPER_MOCK_H_

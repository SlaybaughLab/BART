#ifndef BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_
#define BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_

/*
  deal.II defines a global Assert() macro that clashes with a
  gmock_internal_utils function.
  We want to keep the Assert macro out of
  gmock/internal/gmock-internal-utils.h, but still have it
  available for internal use by any deal.II file that may be
  included after this file.
  
  If Assert is defined, save its value, undefine it, include gmock, and
  redifine it; Else, simply include gmock.
*/

#ifdef Assert
#pragma push_macro("Assert")
#undef Assert
#include "gmock/gmock.h"
#pragma pop_macro("Assert")
#else
#include "gmock/gmock.h"
#endif /*ifdef Assert*/

#endif // BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_

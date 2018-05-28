#ifndef BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_
#define BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_

/*
  deal.II defines a global Assert() macro that clashes with a
  gmock_internal_utils function.
  We want to keep the Assert macro out of
  gmock/internal/gmock-internal-utils.h, but still have it
  available for internal use by any deal.II file that may be
  included after this file.
  
  If Assert is defined, undefine it, include gmock, and
  redifine it; Else, simply include gmock.
*/

#ifdef Assert
#undef Assert
#include "gmock/gmock.h"
/**
	copy/paste of macro from deal.II-v9.0.0/include/deal.II/base/exceptions.h
 */
#ifdef DEBUG
#  ifdef DEAL_II_HAVE_BUILTIN_EXPECT
#    define Assert(cond, exc)                                                \
{                                                                            \
  if (__builtin_expect(!(cond), false))                                      \
    ::dealii::deal_II_exceptions::internals:: issue_error_noreturn(          \
        ::dealii::deal_II_exceptions::internals::abort_on_exception,         \
        __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc);          \
}
#  else /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/
#    define Assert(cond, exc)                                                \
{                                                                            \
  if (!(cond))                                                               \
    ::dealii::deal_II_exceptions::internals:: issue_error_noreturn(          \
        ::dealii::deal_II_exceptions::internals::abort_on_exception,         \
        __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, #exc, exc);          \
}
#  endif /*ifdef DEAL_II_HAVE_BUILTIN_EXPECT*/
#else
#define Assert(cond, exc)                                                    \
  {}
#endif //end of copy/paste
#else
#include "gmock/gmock.h"
#endif /*ifdef Assert*/

#endif // BART_SRC_TEST_HELPERS_GMOCK_WRAPPER_H_

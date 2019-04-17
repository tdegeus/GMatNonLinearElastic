/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

================================================================================================= */

#ifndef GMATNONLINEARELASTIC_H
#define GMATNONLINEARELASTIC_H

// -------------------------------------------------------------------------------------------------

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

// -------------------------------------------------------------------------------------------------

#ifndef NDEBUG
#define GMATNONLINEARELASTIC_ENABLE_ASSERT
#endif

#ifdef GMATNONLINEARELASTIC_ENABLE_ASSERT
#define GMATNONLINEARELASTIC_ASSERT(expr) GMATNONLINEARELASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
#define GMATNONLINEARELASTIC_ASSERT_IMPL(expr, file, line)                                                                \
    if (!(expr))                                                                                                          \
    {                                                                                                                     \
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": assertion failed (" #expr ") \n\t"); \
    }
#else
#define GMATNONLINEARELASTIC_ASSERT(expr)
#endif

// -------------------------------------------------------------------------------------------------

#define GMATNONLINEARELASTIC_WORLD_VERSION 0
#define GMATNONLINEARELASTIC_MAJOR_VERSION 0
#define GMATNONLINEARELASTIC_MINOR_VERSION 3

#define GMATNONLINEARELASTIC_VERSION_AT_LEAST(x,y,z) \
  (GMATNONLINEARELASTIC_WORLD_VERSION>x || (GMATNONLINEARELASTIC_WORLD_VERSION>=x && \
  (GMATNONLINEARELASTIC_MAJOR_VERSION>y || (GMATNONLINEARELASTIC_MAJOR_VERSION>=y && \
                                            GMATNONLINEARELASTIC_MINOR_VERSION>=z))))

#define GMATNONLINEARELASTIC_VERSION(x,y,z) \
  (GMATNONLINEARELASTIC_WORLD_VERSION==x && \
   GMATNONLINEARELASTIC_MAJOR_VERSION==y && \
   GMATNONLINEARELASTIC_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#endif

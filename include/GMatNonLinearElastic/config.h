/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_H
#define GMATNONLINEARELASTIC_H

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

#ifdef GMATNONLINEARELASTIC_ENABLE_ASSERT

    #define GMATNONLINEARELASTIC_ASSERT(expr) \
        GMATNONLINEARELASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)

    #define GMATNONLINEARELASTIC_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }

#else

    #define GMATNONLINEARELASTIC_ASSERT(expr)

#endif

#define GMATNONLINEARELASTIC_VERSION_MAJOR 0
#define GMATNONLINEARELASTIC_VERSION_MINOR 1
#define GMATNONLINEARELASTIC_VERSION_PATCH 0

#define GMATNONLINEARELASTIC_VERSION_AT_LEAST(x,y,z) \
    (GMATNONLINEARELASTIC_VERSION_MAJOR > x || (GMATNONLINEARELASTIC_VERSION_MAJOR >= x && \
    (GMATNONLINEARELASTIC_VERSION_MINOR > y || (GMATNONLINEARELASTIC_VERSION_MINOR >= y && \
                                                GMATNONLINEARELASTIC_VERSION_PATCH >= z))))

#define GMATNONLINEARELASTIC_VERSION(x,y,z) \
    (GMATNONLINEARELASTIC_VERSION_MAJOR == x && \
     GMATNONLINEARELASTIC_VERSION_MINOR == y && \
     GMATNONLINEARELASTIC_VERSION_PATCH == z)

#endif

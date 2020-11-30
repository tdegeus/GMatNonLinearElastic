/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_H
#define GMATNONLINEARELASTIC_H

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

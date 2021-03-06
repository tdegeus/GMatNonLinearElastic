/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_HPP
#define GMATNONLINEARELASTIC_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatNonLinearElastic {
namespace Cartesian3d {

template <class T, class U>
inline void epseq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(2.0 / 3.0);
}

template <class T>
inline auto Epseq(const T& A)
{
    return xt::eval(std::sqrt(2.0 / 3.0) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

template <class T, class U>
inline void sigeq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(1.5);
}

template <class T>
inline auto Sigeq(const T& A)
{
    return xt::eval(std::sqrt(1.5) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

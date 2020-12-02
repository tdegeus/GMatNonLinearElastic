/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_H
#define GMATNONLINEARELASTIC_CARTESIAN3D_H

#include <GMatTensor/Cartesian3d.h>

#include "config.h"

namespace GMatNonLinearElastic {
namespace Cartesian3d {

// Unit tensors

using GMatTensor::Cartesian3d::I2;
using GMatTensor::Cartesian3d::II;
using GMatTensor::Cartesian3d::I4;
using GMatTensor::Cartesian3d::I4rt;
using GMatTensor::Cartesian3d::I4s;
using GMatTensor::Cartesian3d::I4d;

// Tensor decomposition

using GMatTensor::Cartesian3d::hydrostatic;
using GMatTensor::Cartesian3d::Hydrostatic;
using GMatTensor::Cartesian3d::deviatoric;
using GMatTensor::Cartesian3d::Deviatoric;

// Equivalent strain

template <class T, class U>
inline void epseq(const T& A, U& ret);

template <class T>
inline auto Epseq(const T& A);

// Equivalent stress

template <class T, class U>
inline void sigeq(const T& A, U& ret);

template <class T>
inline auto Sigeq(const T& A);

// Material point

class NonLinearElastic
{
public:
    NonLinearElastic() = default;
    NonLinearElastic(double kappa, double sig0, double eps0, double m);

    double kappa() const;
    double sig0() const;
    double eps0() const;
    double m() const;

    template <class T> void setStrain(const T& arg, bool compute_tangent = true);
    template <class T> void strain(T& ret) const;
    template <class T> void stress(T& ret) const;
    template <class T> void tangent(T& ret) const;

    template <class T> void setStrainIterator(const T* arg, bool compute_tangent = true);
    template <class T> void strainIterator(T* ret) const;
    template <class T> void stressIterator(T* ret) const;
    template <class T> void tangentIterator(T* ret) const;

    xt::xtensor<double, 2> Strain() const;
    xt::xtensor<double, 2> Stress() const;
    xt::xtensor<double, 4> Tangent() const;

private:
    double m_kappa;
    double m_sig0;
    double m_eps0;
    double m_m;
    std::array<double, 9> m_Eps; // strain tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_Sig; // stress tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    xt::xtensor<double, 4> m_C; // tangent stiffness
};

// Material identifier

struct Type {
    enum Value {
        Unset,
        NonLinearElastic,
    };
};

// Array of material points

template <size_t N>
class Array : public GMatTensor::Cartesian3d::Array<N>
{
public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    // Constructors

    Array() = default;
    Array(const std::array<size_t, N>& shape);
    Array(const std::array<size_t, N>& shape, double kappa, double sig0, double eps0, double m);

    // Overloaded methods:
    // - "shape"
    // - unit tensors: "I2", "II", "I4", "I4rt", "I4s", "I4d"

    // Type

    xt::xtensor<size_t, N> type() const;

    // Parameters

    xt::xtensor<double, N> kappa() const;
    xt::xtensor<double, N> sig0() const;
    xt::xtensor<double, N> eps0() const;
    xt::xtensor<double, N> m() const;

    // Set parameters for a batch of points

    void setNonLinearElastic(
        const xt::xtensor<size_t, N>& I,
        double kappa,
        double sig0,
        double eps0,
        double m);

    // Set strain tensor, get the response

    void setStrain(const xt::xtensor<double, N + 2>& arg, bool compute_tangent = true);
    void strain(xt::xtensor<double, N + 2>& ret) const;
    void stress(xt::xtensor<double, N + 2>& ret) const;
    void tangent(xt::xtensor<double, N + 4>& ret) const;

    // Auto-allocation of the functions above

    xt::xtensor<double, N + 2> Strain() const;
    xt::xtensor<double, N + 2> Stress() const;
    xt::xtensor<double, N + 4> Tangent() const;

    // Get copy or reference to the underlying model at on point

    auto getNonLinearElastic(const std::array<size_t, N>& index) const;
    auto* refNonLinearElastic(const std::array<size_t, N>& index);

private:
    // Material vectors
    std::vector<NonLinearElastic> m_NonLinearElastic;

    // Identifiers for each matrix entry
    xt::xtensor<size_t, N> m_type;  // type (e.g. "Type::Elastic")
    xt::xtensor<size_t, N> m_index; // index from the relevant material vector (e.g. "m_Elastic")

    // Shape
    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;
};

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#include "Cartesian3d.hpp"
#include "Cartesian3d_Array.hpp"
#include "Cartesian3d_NonLinearElastic.hpp"

#endif

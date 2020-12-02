/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_MATRIX_HPP
#define GMATNONLINEARELASTIC_CARTESIAN3D_MATRIX_HPP

#include "Cartesian3d.h"

namespace GMatNonLinearElastic {
namespace Cartesian3d {

template <size_t N>
inline Array<N>::Array(const std::array<size_t, N>& shape)
{
    this->init(shape);
    m_type = xt::ones<size_t>(m_shape) * Type::Unset;
    m_index = xt::empty<size_t>(m_shape);
}

template <size_t N>
inline Array<N>::Array(
    const std::array<size_t, N>& shape,
    double kappa,
    double sig0,
    double eps0,
    double m)
{
    this->init(shape);
    m_type = xt::ones<size_t>(m_shape) * Type::NonLinearElastic;
    m_index = xt::arange<size_t>(m_size).reshape(m_shape);

    for (size_t i = 0; i < m_size; ++i) {
        m_NonLinearElastic.push_back(NonLinearElastic(kappa, sig0, eps0, m));
    }
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::kappa() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::NonLinearElastic:
            ret.data()[i] = m_NonLinearElastic[m_index.data()[i]].kappa();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::sig0() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::NonLinearElastic:
            ret.data()[i] = m_NonLinearElastic[m_index.data()[i]].sig0();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::eps0() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::NonLinearElastic:
            ret.data()[i] = m_NonLinearElastic[m_index.data()[i]].eps0();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<double, N> Array<N>::m() const
{
    xt::xtensor<double, N> ret = xt::empty<double>(m_shape);

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            ret.data()[i] = 0.0;
            break;
        case Type::NonLinearElastic:
            ret.data()[i] = m_NonLinearElastic[m_index.data()[i]].m();
            break;
        }
    }

    return ret;
}

template <size_t N>
inline xt::xtensor<size_t, N> Array<N>::type() const
{
    return m_type;
}

template <size_t N>
inline void Array<N>::setNonLinearElastic(
    const xt::xtensor<size_t, N>& I,
    double kappa,
    double sig0,
    double eps0,
    double m)
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(m_type, I.shape()));
    GMATNONLINEARELASTIC_ASSERT(xt::all(xt::equal(I, 0ul) || xt::equal(I, 1ul)));
    GMATNONLINEARELASTIC_ASSERT(
        xt::all(xt::equal(xt::where(xt::equal(I, 1ul), m_type, Type::Unset), Type::Unset)));

    for (size_t i = 0; i < m_size; ++i) {
        if (I.data()[i] == 1ul) {
            m_type.data()[i] = Type::NonLinearElastic;
            m_index.data()[i] = m_NonLinearElastic.size();
            m_NonLinearElastic.push_back(NonLinearElastic(kappa, sig0, eps0, m));
        }
    }
}

template <size_t N>
inline void Array<N>::setStrain(const xt::xtensor<double, N + 2>& A, bool tangent)
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            break;
        case Type::NonLinearElastic:
            m_NonLinearElastic[m_index.data()[i]].setStrainIterator(
                &A.data()[i * m_stride_tensor2],
                tangent);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::strain(xt::xtensor<double, N + 2>& A) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O2(&A.data()[i * m_stride_tensor2]);
            break;
        case Type::NonLinearElastic:
            m_NonLinearElastic[m_index.data()[i]].strainIterator(&A.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::stress(xt::xtensor<double, N + 2>& A) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor2));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O2(&A.data()[i * m_stride_tensor2]);
            break;
        case Type::NonLinearElastic:
            m_NonLinearElastic[m_index.data()[i]].stressIterator(&A.data()[i * m_stride_tensor2]);
            break;
        }
    }
}

template <size_t N>
inline void Array<N>::tangent(xt::xtensor<double, N + 4>& A) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(A, m_shape_tensor4));

    #pragma omp parallel for
    for (size_t i = 0; i < m_size; ++i) {
        switch (m_type.data()[i]) {
        case Type::Unset:
            GMatTensor::Cartesian3d::pointer::O4(&A.data()[i * m_stride_tensor4]);
            break;
        case Type::NonLinearElastic:
            m_NonLinearElastic[m_index.data()[i]].tangentIterator(&A.data()[i * m_stride_tensor4]);
            break;
        }
    }
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Strain() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->strain(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 2> Array<N>::Stress() const
{
    xt::xtensor<double, N + 2> ret = xt::empty<double>(m_shape_tensor2);
    this->stress(ret);
    return ret;
}

template <size_t N>
inline xt::xtensor<double, N + 4> Array<N>::Tangent() const
{
    xt::xtensor<double, N + 4> ret = xt::empty<double>(m_shape_tensor4);
    this->tangent(ret);
    return ret;
}

template <size_t N>
inline auto Array<N>::getNonLinearElastic(const std::array<size_t, N>& index) const
{
    GMATNONLINEARELASTIC_ASSERT(m_type[index] == Type::NonLinearElastic);
    return m_NonLinearElastic[m_index[index]];
}

template <size_t N>
inline auto* Array<N>::refNonLinearElastic(const std::array<size_t, N>& index)
{
    GMATNONLINEARELASTIC_ASSERT(m_type[index] == Type::NonLinearElastic);
    return &m_NonLinearElastic[m_index[index]];
}

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

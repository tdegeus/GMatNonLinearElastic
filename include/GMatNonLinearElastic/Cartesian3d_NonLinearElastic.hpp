/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

*/

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_NONLINEARELASTIC_HPP
#define GMATNONLINEARELASTIC_CARTESIAN3D_NONLINEARELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatNonLinearElastic {
namespace Cartesian3d {

inline NonLinearElastic::NonLinearElastic(double kappa, double sig0, double eps0, double m)
    : m_kappa(kappa), m_sig0(sig0), m_eps0(eps0), m_m(m)
{
    m_C = xt::empty<double>({3, 3, 3, 3});
}

inline double NonLinearElastic::kappa() const
{
    return m_kappa;
}

inline double NonLinearElastic::sig0() const
{
    return m_sig0;
}

inline double NonLinearElastic::eps0() const
{
    return m_eps0;
}

inline double NonLinearElastic::m() const
{
    return m_m;
}

template <class T>
inline void NonLinearElastic::setStrainIterator(const T* arg, bool tangent)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;
    std::copy(arg, arg + 9, m_Eps.begin());

    std::array<double, 9> Epsd;
    double epsm = GT::hydrostatic_deviatoric(&m_Eps[0], &Epsd[0]);
    double epseq = std::sqrt(2.0 / 3.0 * GT::A2s_ddot_B2s(&Epsd[0], &Epsd[0]));

    if (epseq != 0.0) {
        double f = 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m) * std::pow(epseq, m_m - 1.0);
        m_Sig[0] = f * Epsd[0] + 3.0 * m_kappa * epsm;
        m_Sig[1] = f * Epsd[1];
        m_Sig[2] = f * Epsd[2];
        m_Sig[3] = f * Epsd[3];
        m_Sig[4] = f * Epsd[4] + 3.0 * m_kappa * epsm;
        m_Sig[5] = f * Epsd[5];
        m_Sig[6] = f * Epsd[6];
        m_Sig[7] = f * Epsd[7];
        m_Sig[8] = f * Epsd[8] + 3.0 * m_kappa * epsm;
    }

    m_Sig[0] = 3.0 * m_kappa * epsm;
    m_Sig[1] = 0.0;
    m_Sig[2] = 0.0;
    m_Sig[3] = 0.0;
    m_Sig[4] = 3.0 * m_kappa * epsm;
    m_Sig[5] = 0.0;
    m_Sig[6] = 0.0;
    m_Sig[7] = 0.0;
    m_Sig[8] = 3.0 * m_kappa * epsm;

    if (!tangent) {
        return;
    }

    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();

    if (epseq != 0.0) {
        GMatTensor::Cartesian3d::pointer::A2_dyadic_B2(&Epsd[0], &Epsd[0], m_C.data());
        m_C *= 2.0 / 3.0 * (m_m - 1.0) * std::pow(epseq, m_m - 3.0);
        m_C += std::pow(epseq, m_m - 1.0) * I4d;
        m_C *= 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m);
        m_C += m_kappa * II;
    }
    else {
        xt::noalias(m_C) = m_kappa * II + I4d;
    }
}

template <class T>
inline void NonLinearElastic::strainIterator(T* ret) const
{
    std::copy(m_Eps.begin(), m_Eps.end(), ret);
}

template <class T>
inline void NonLinearElastic::stressIterator(T* ret) const
{
    std::copy(m_Sig.begin(), m_Sig.end(), ret);
}

template <class T>
inline void NonLinearElastic::tangentIterator(T* ret) const
{
    std::copy(m_C.cbegin(), m_C.cend(), ret);
}

template <class T>
inline void NonLinearElastic::setStrain(const T& a, bool tangent)
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->setStrainIterator(a.data(), tangent);
}

template <class T>
inline void NonLinearElastic::strain(T& a) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->strainIterator(a.data());
}

template <class T>
inline void NonLinearElastic::stress(T& a) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->stressIterator(a.data());
}

template <class T>
inline void NonLinearElastic::tangent(T& a) const
{
    GMATNONLINEARELASTIC_ASSERT(xt::has_shape(a, {3, 3, 3, 3}));
    return this->stressIterator(a.data());
}

inline xt::xtensor<double, 2> NonLinearElastic::Strain() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->strainIterator(ret.data());
    return ret;
}

inline xt::xtensor<double, 2> NonLinearElastic::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->stressIterator(ret.data());
    return ret;
}

inline xt::xtensor<double, 4> NonLinearElastic::Tangent() const
{
    xt::xtensor<double, 4> ret = xt::empty<double>({3, 3, 3, 3});
    this->tangentIterator(ret.data());
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

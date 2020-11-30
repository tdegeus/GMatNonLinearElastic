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
inline void NonLinearElastic::setStrain(const T& a)
{
    GMATELASTIC_ASSERT(xt::has_shape(a, {3, 3}));
    return this->setStrainIterator(a.cbegin());
}

template <class T>
inline void NonLinearElastic::setStrainIterator(const T& begin)
{
    std::copy(begin, begin + 9, m_Eps.begin());

    double epsm = GMatTensor::Cartesian3d::pointer::hydrostatic_deviatoric(&m_Eps[0], &m_Epsd[0]);
    m_epseq = std::sqrt(2.0 / 3.0 * GMatTensor::Cartesian3d::pointer::A2_ddot_B2(&m_Epsd[0], &m_Epsd[0]));

    if (m_epseq != 0.0) {
        double f = 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m) * std::pow(m_epseq, m_m - 1.0);
        m_Sig[0] = f * m_Epsd[0] + 3.0 * m_kappa * epsm;
        m_Sig[1] = f * m_Epsd[1];
        m_Sig[2] = f * m_Epsd[2];
        m_Sig[3] = f * m_Epsd[3];
        m_Sig[4] = f * m_Epsd[4] + 3.0 * m_kappa * epsm;
        m_Sig[5] = f * m_Epsd[5];
        m_Sig[6] = f * m_Epsd[6];
        m_Sig[7] = f * m_Epsd[7];
        m_Sig[8] = f * m_Epsd[8] + 3.0 * m_kappa * epsm;
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
}

template <class T>
inline void NonLinearElastic::strain(T& a) const
{
    GMATELASTOPLASTICQPOT_ASSERT(xt::has_shape(a, {2, 2}));
    return this->strainIterator(a.begin());
}

template <class T>
inline void NonLinearElastic::strainIterator(const T& begin) const
{
    std::copy(m_Eps.begin(), m_Eps.end(), begin);
}

inline xt::xtensor<double, 2> NonLinearElastic::Strain() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({2, 2});
    this->strainIterator(ret.begin());
    return ret;
}

template <class T>
inline void NonLinearElastic::stress(T& a) const
{
    GMATELASTOPLASTICQPOT_ASSERT(xt::has_shape(a, {2, 2}));
    return this->stressIterator(a.begin());
}

template <class T>
inline void NonLinearElastic::stressIterator(const T& begin) const
{
    std::copy(m_Sig.begin(), m_Sig.end(), begin);
}

inline xt::xtensor<double, 2> NonLinearElastic::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({2, 2});
    this->stressIterator(ret.begin());
    return ret;
}

template <class T>
inline void NonLinearElastic::tangent(T& C) const
{
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();

    if (m_epseq != 0.0) {
        xt::xtensor<double, 2> K = xt::empty<double>({3, 3, 3, 3});
        xt::xtensor<double, 2> Epsd = {
            {m_Epsd[0], m_Epsd[1], m_Epsd[2]},
            {m_Epsd[3], m_Epsd[4], m_Epsd[5]},
            {m_Epsd[6], m_Epsd[7], m_Epsd[8]}};
        GMatTensor::Cartesian3d::xtensor::A2_dyadic_B2(Epsd, Epsd, K);
        K *= 2.0 / 3.0 * (m_m - 1.0) * std::pow(m_epseq, m_m - 3.0);
        K += std::pow(m_epseq, m_m - 1.0) * I4d;
        K *= 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m);
        K += m_kappa * II;
        xt::noalias(C) = K;
    }
    else {
        xt::noalias(C) = m_kappa * II + I4d;
    }
}

inline xt::xtensor<double, 4> NonLinearElastic::Tangent() const
{
    xt::xtensor<double, 4> ret = xt::zeros<double>({3, 3, 3, 3});
    this->tangent(ret);
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

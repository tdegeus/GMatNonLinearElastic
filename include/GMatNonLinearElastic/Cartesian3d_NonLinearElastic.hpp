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
inline void NonLinearElastic::stress(const Tensor2& Eps, T&& Sig) const
{
    auto I = Cartesian3d::I2();
    auto epsm = trace(Eps) / 3.0;
    auto Epsd = Eps - epsm * I;
    auto epseq = std::sqrt(2.0 / 3.0 * A2_ddot_B2(Epsd, Epsd));
    ;

    if (epseq != 0.0) {
        xt::noalias(Sig)
            = 3.0 * m_kappa * epsm * I
            + 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m) * std::pow(epseq, m_m - 1.0) * Epsd;
    }
    else {
        xt::noalias(Sig) = 3.0 * m_kappa * epsm * I;
    }
}

inline Tensor2 NonLinearElastic::Stress(const Tensor2& Eps) const
{
    Tensor2 Sig;
    this->stress(Eps, Sig);
    return Sig;
}

template <class T, class S>
inline void NonLinearElastic::tangent(const Tensor2& Eps, T&& Sig, S&& C) const
{
    auto I = Cartesian3d::I2();
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto epsm = trace(Eps) / 3.0;
    auto Epsd = Eps - epsm * I;
    auto epseq = std::sqrt(2.0 / 3.0 * A2_ddot_B2(Epsd, Epsd));

    if (epseq != 0.0) {
        xt::noalias(Sig)
            = 3.0 * m_kappa * epsm * I
            + 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m) * std::pow(epseq, m_m - 1.0) * Epsd;
    }
    else {
        xt::noalias(Sig) = 3.0 * m_kappa * epsm * I;
    }

    if (epseq != 0.0) {
        Tensor4 K;
        A2_dyadic_B2(Epsd, Epsd, K);
        K *= 2.0 / 3.0 * (m_m - 1.0) * std::pow(epseq, m_m - 3.0);
        K += std::pow(epseq, m_m - 1.0) * I4d;
        K *= 2.0 / 3.0 * m_sig0 / std::pow(m_eps0, m_m);
        K += m_kappa * II;
        xt::noalias(C) = K;
    }
    else {
        xt::noalias(C) = m_kappa * II + I4d;
    }
}

inline std::tuple<Tensor2, Tensor4> NonLinearElastic::Tangent(const Tensor2& Eps) const
{
    Tensor2 Sig;
    Tensor4 C;
    this->tangent(Eps, Sig, C);
    return std::make_tuple(Sig, C);
}

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

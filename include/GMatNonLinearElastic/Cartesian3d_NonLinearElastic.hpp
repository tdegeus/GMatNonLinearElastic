/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

================================================================================================= */

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_NONLINEARELASTIC_HPP
#define GMATNONLINEARELASTIC_CARTESIAN3D_NONLINEARELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatNonLinearElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline NonLinearElastic::NonLinearElastic(double kappa, double sig0, double eps0, double m) :
  m_kappa(kappa), m_sig0(sig0), m_eps0(eps0), m_m(m)
{
}

// -------------------------------------------------------------------------------------------------

inline double NonLinearElastic::kappa() const
{
  return m_kappa;
}

// -------------------------------------------------------------------------------------------------

inline double NonLinearElastic::sig0() const
{
  return m_sig0;
}

// -------------------------------------------------------------------------------------------------

inline double NonLinearElastic::eps0() const
{
  return m_eps0;
}

// -------------------------------------------------------------------------------------------------

inline double NonLinearElastic::m() const
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void NonLinearElastic::stress(const T2& Eps, T&& Sig) const
{
  auto I     = Cartesian3d::I();
  auto epsm  = trace(Eps) / 3.0;
  auto Epsd  = Eps - epsm * I;
  auto epseq = std::sqrt(2.0/3.0 * ddot22(Epsd,Epsd));;

  if ( epseq != 0.0 )
  {
    xt::noalias(Sig) = 3.*m_kappa*epsm*I + 2./3.*m_sig0/std::pow(m_eps0,m_m)*std::pow(epseq,m_m-1.)*Epsd;
  }
  else
  {
    xt::noalias(Sig) = 3.*m_kappa*epsm*I;
  }
}

// -------------------------------------------------------------------------------------------------

inline T2 NonLinearElastic::Stress(const T2& Eps) const
{
  T2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void NonLinearElastic::tangent(const T2& Eps, T&& Sig, S&& C) const
{
  auto I     = Cartesian3d::I();
  auto II    = Cartesian3d::II();
  auto I4d   = Cartesian3d::I4d();
  auto epsm  = trace(Eps) / 3.0;
  auto Epsd  = Eps - epsm * I;
  auto epseq = std::sqrt(2.0/3.0 * ddot22(Epsd,Epsd));

  if ( epseq != 0.0 )
  {
    xt::noalias(Sig) = 3.*m_kappa*epsm*I + 2./3.*m_sig0/std::pow(m_eps0,m_m)*std::pow(epseq,m_m-1.)*Epsd;
  }
  else
  {
    xt::noalias(Sig) = 3.*m_kappa*epsm*I;
  }

  if ( epseq != 0.0 )
  {
    T4 out;
    dyadic22(Epsd, Epsd, out);
    out *= 2./3. * (m_m-1.) * std::pow(epseq,m_m-3.);
    out += std::pow(epseq, m_m-1.) * I4d;
    out *= 2./3. * m_sig0 / std::pow(m_eps0,m_m);
    out += m_kappa * II;
    xt::noalias(C) = out;
  }
  else
  {
    xt::noalias(C) = m_kappa * II + I4d;
  }
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> NonLinearElastic::Tangent(const T2& Eps) const
{
  T2 Sig;
  T4 C;
  this->tangent(Eps, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif

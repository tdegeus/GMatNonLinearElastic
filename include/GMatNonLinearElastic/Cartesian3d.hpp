/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot

================================================================================================= */

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_HPP
#define GMATNONLINEARELASTIC_CARTESIAN3D_HPP

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.h"

// =================================================================================================

namespace GMatNonLinearElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

template<class T>
inline double trace(const T &A)
{
  return A(0,0) + A(1,1) + A(2,2);
}

// -------------------------------------------------------------------------------------------------

template<class T>
inline double ddot22(const T &A, const T &B)
{
  return A(0,0) * B(0,0) + 2.0 * A(0,1) * B(0,1) + 2.0 * A(0,2) * B(0,2) +
         A(1,1) * B(1,1) + 2.0 * A(1,2) * B(1,2) +
         A(2,2) * B(2,2);
}

// -------------------------------------------------------------------------------------------------

template<class T2, class T4>
inline void dyadic22(const T2 &A, const T2 &B, T4 &C)
{
  C.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);
}

// -------------------------------------------------------------------------------------------------

inline T2 I()
{
  return T2({{1., 0., 0.},
             {0., 1., 0.},
             {0., 0., 1.}});
}

// -------------------------------------------------------------------------------------------------

inline T4 II()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == j and k == l )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == l and j == k )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4rt()
{
  T4 out;

  out.fill(0.0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          if ( i == k and j == l )
            out(i,j,k,l) = 1.;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline T4 I4s()
{
  return .5 * ( I4() + I4rt() );
}

// -------------------------------------------------------------------------------------------------

inline T4 I4d()
{
  return I4s() - II()/3.;
}

// -------------------------------------------------------------------------------------------------

inline T2 NonLinearElastic::Sig(const T2 &Eps) const
{
  // define identity tensor
  T2 I = Cartesian3d::I();

  // decompose strain, compute equivalent strain
  double epsm  = trace(Eps)/3.;
  auto   Epsd  = Eps - epsm * I;
  double epseq = std::sqrt(2./3.*ddot22(Epsd,Epsd));;

  // compute stress, guard for zero-division
  if ( epseq != 0.0 )
  {
    return 3.*m_kappa*epsm*I + 2./3.*m_sig0/std::pow(m_eps0,m_m)*std::pow(epseq,m_m-1.)*Epsd;
  }
  else
  {
    return 3.*m_kappa*epsm*I;
  }
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<T2,T4> NonLinearElastic::Tangent(const T2 &Eps) const
{
  // define identity tensor, allocate result
  T2 I   = Cartesian3d::I();
  T4 II  = Cartesian3d::II();
  T4 I4d = Cartesian3d::I4d();
  T2 Sig;
  T4 C4;

  // decompose strain, compute equivalent strain
  double epsm  = trace(Eps)/3.;
  auto   Epsd  = Eps - epsm * I;
  double epseq = std::sqrt(2./3.*ddot22(Epsd,Epsd));;

  // compute stress, guard for zero-division
  if ( epseq != 0.0 )
  {
    Sig = 3.*m_kappa*epsm*I + 2./3.*m_sig0/std::pow(m_eps0,m_m)*std::pow(epseq,m_m-1.)*Epsd;
  }
  else
  {
    Sig = 3.*m_kappa*epsm*I;
  }

  // compute tangent, guard for zero-division
  if ( epseq != 0.0 )
  {
    dyadic22(Epsd,Epsd,C4);
    C4 *= 2./3. * (m_m-1.) * std::pow(epseq,m_m-3.);
    C4 += std::pow(epseq,m_m-1.) * I4d;
    C4 *= 2./3. * m_sig0 / std::pow(m_eps0,m_m);
    C4 += m_kappa * II;
  }
  else
  {
    C4 = m_kappa * II + I4d;
  }

  return std::make_tuple(Sig, C4);
}

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip) : m_nelem(nelem), m_nip(nip)
{
  m_set   = xt::zeros<int   >({nelem, nip});
  m_kappa = xt::empty<double>({nelem, nip});
  m_sig0  = xt::empty<double>({nelem, nip});
  m_eps0  = xt::empty<double>({nelem, nip});
  m_m     = xt::empty<double>({nelem, nip});
}

// -------------------------------------------------------------------------------------------------

inline Matrix::Matrix(size_t nelem, size_t nip, double kappa, double sig0, double eps0, double m) :
  m_nelem(nelem), m_nip(nip)
{
  m_set   = xt::ones<int   >({nelem, nip});
  m_kappa = xt::ones<double>({nelem, nip}) * kappa;
  m_sig0  = xt::ones<double>({nelem, nip}) * sig0;
  m_eps0  = xt::ones<double>({nelem, nip}) * eps0;
  m_m     = xt::ones<double>({nelem, nip}) * m;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::check() const
{
  for ( size_t e = 0 ; e < m_nelem ; ++e )
    for ( size_t q = 0 ; q < m_nip ; ++q )
      if ( m_set(e,q) == 0 )
        throw std::runtime_error("No type set for: "+std::to_string(e)+", "+std::to_string(q));
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::set(
  const xt::xtensor<size_t,2> &I, double kappa, double sig0, double eps0, double m)
{
  // check input
  #ifndef NDEBUG
    // - shape
    assert( I.shape() == m_set.shape() );
    // - empty material definitions
    for ( size_t e = 0 ; e < m_nelem ; ++e )
      for ( size_t q = 0 ; q < m_nip ; ++q )
        if ( I(e,q) )
          assert( m_set(e,q) == 0 );
  #endif

  // set type and position in material vector
  for ( size_t e = 0 ; e < m_nelem ; ++e ) {
    for ( size_t q = 0 ; q < m_nip ; ++q ) {
      if ( I(e,q) ) {
        m_set  (e,q) = 1;
        m_kappa(e,q) = kappa;
        m_sig0 (e,q) = sig0;
        m_eps0 (e,q) = eps0;
        m_m    (e,q) = m;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::I() const
{
  xt::xtensor<double,4> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T2 unit = Cartesian3d::I();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0), xt::xshape<m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::II() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::II();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4rt() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4rt();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4s() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4s();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,6> Matrix::I4d() const
{
  xt::xtensor<double,6> out = xt::empty<double>({m_nelem, m_nip, m_ndim, m_ndim, m_ndim, m_ndim});

  #pragma omp parallel
  {
    T4 unit = Cartesian3d::I4d();

    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        auto view = xt::adapt(&out(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());

        xt::noalias(view) = unit;
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::Sig(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Sig) const
{
  assert( a_Eps.shape()[0] == m_nelem       );
  assert( a_Eps.shape()[1] == m_nip         );
  assert( a_Eps.shape()[2] == m_ndim        );
  assert( a_Eps.shape()[3] == m_ndim        );
  assert( a_Eps.shape()    == a_Sig.shape() );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensor
    T2 I = Cartesian3d::I();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        // alias
        auto  Eps   = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig   = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto& kappa = m_kappa(e,q);
        auto& sig0  = m_sig0 (e,q);
        auto& eps0  = m_eps0 (e,q);
        auto& m     = m_m    (e,q);

        // decompose strain
        double epsm  = trace(Eps)/3.;
        auto   Epsd  = Eps - epsm * I;
        double epseq = std::sqrt(2./3.*ddot22(Epsd,Epsd));;

        // compute stress
        if ( epseq != 0.0 )
          xt::noalias(Sig) = 3.*kappa*epsm*I + 2./3.*sig0/std::pow(eps0,m)*std::pow(epseq,m-1.)*Epsd;
        else
          xt::noalias(Sig) = 3.*kappa*epsm*I;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void Matrix::Tangent(const xt::xtensor<double,4> &a_Eps,
  xt::xtensor<double,4> &a_Sig, xt::xtensor<double,6> &a_Tangent) const
{
  assert( a_Eps.shape()[0]     == m_nelem       );
  assert( a_Eps.shape()[1]     == m_nip         );
  assert( a_Eps.shape()[2]     == m_ndim        );
  assert( a_Eps.shape()[3]     == m_ndim        );
  assert( a_Eps.shape()        == a_Sig.shape() );
  assert( a_Tangent.shape()[0] == m_nelem       );
  assert( a_Tangent.shape()[1] == m_nip         );
  assert( a_Tangent.shape()[2] == m_ndim        );
  assert( a_Tangent.shape()[3] == m_ndim        );
  assert( a_Tangent.shape()[4] == m_ndim        );
  assert( a_Tangent.shape()[5] == m_ndim        );

  // start threads (all allocated variables inside this block are local to each thread)
  #pragma omp parallel
  {
    // identity tensors
    T2 I   = Cartesian3d::I();
    T4 II  = Cartesian3d::II();
    T4 I4d = Cartesian3d::I4d();

    // loop over all points
    #pragma omp for
    for ( size_t e = 0 ; e < m_nelem ; ++e )
    {
      for ( size_t q = 0 ; q < m_nip ; ++q )
      {
        // alias
        auto  Eps   = xt::adapt(&a_Eps(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  Sig   = xt::adapt(&a_Sig(e,q,0,0), xt::xshape<m_ndim,m_ndim>());
        auto  C4    = xt::adapt(&a_Tangent(e,q,0,0,0,0), xt::xshape<m_ndim,m_ndim,m_ndim,m_ndim>());
        auto& kappa = m_kappa(e,q);
        auto& sig0  = m_sig0 (e,q);
        auto& eps0  = m_eps0 (e,q);
        auto& m     = m_m    (e,q);

        // decompose strain
        double epsm  = trace(Eps)/3.;
        auto   Epsd  = Eps - epsm * I;
        double epseq = std::sqrt(2./3.*ddot22(Epsd,Epsd));;

        // compute stress
        if ( epseq != 0.0 )
        {
          xt::noalias(Sig) = 3.*kappa*epsm*I + 2./3.*sig0/std::pow(eps0,m)*std::pow(epseq,m-1.)*Epsd;
        }
        else
        {
          xt::noalias(Sig) = 3.*kappa*epsm*I;
        }

        // compute tangent
        if ( epseq != 0.0 )
        {
          dyadic22(Epsd,Epsd,C4);
          xt::noalias(C4) = C4 * 2./3. * (m-1.) * std::pow(epseq,m-3.);
          xt::noalias(C4) = C4 + std::pow(epseq,m-1.) * I4d;
          xt::noalias(C4) = C4 * 2./3. * sig0 / std::pow(eps0,m);
          xt::noalias(C4) = C4 + kappa * II;
        }
        else
        {
          xt::noalias(C4) = kappa * II + I4d;
        }
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Matrix::Sig(const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,4> a_Sig = xt::empty<double>(a_Eps.shape());

  this->Sig(a_Eps, a_Sig);

  return a_Sig;
}

// -------------------------------------------------------------------------------------------------

inline std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Matrix::Tangent(
  const xt::xtensor<double,4> &a_Eps) const
{
  xt::xtensor<double,4> a_Sig     = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim});
  xt::xtensor<double,6> a_Tangent = xt::empty<double>({m_nelem,m_nip,m_ndim,m_ndim,m_ndim,m_ndim});

  this->Tangent(a_Eps, a_Sig, a_Tangent);

  return std::make_tuple(a_Sig, a_Tangent);
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif

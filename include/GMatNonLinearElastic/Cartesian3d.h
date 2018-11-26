/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GMatNonLinearElastic

================================================================================================= */

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_H
#define GMATNONLINEARELASTIC_CARTESIAN3D_H

// -------------------------------------------------------------------------------------------------

#include "config.h"

// =================================================================================================

namespace GMatNonLinearElastic {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

using T2 = xt::xtensor_fixed<double, xt::xshape<3,3>>;
using T4 = xt::xtensor_fixed<double, xt::xshape<3,3,3,3>>;

template<class T> double trace (const T &A);
template<class T> double ddot22(const T &A, const T &B);

T2 I();
T4 II();
T4 I4();
T4 I4rt();
T4 I4s();
T4 I4d();

// -------------------------------------------------------------------------------------------------

class NonLinearElastic
{
public:

  // constructors
  NonLinearElastic() = default;
  NonLinearElastic(double kappa, double sig0, double eps0, double m);

  // return parameters
  double kappa() const;
  double sig0() const;
  double eps0() const;
  double m() const;

  // compute stress
  T2 Sig(const T2 &Eps) const;

  // compute tangent
  std::tuple<T2,T4> Tangent(const T2 &Eps) const;

private:

  double m_kappa;
  double m_sig0;
  double m_eps0;
  double m_m;
};

// -------------------------------------------------------------------------------------------------

class Matrix
{
public:

  // constructors
  Matrix() = default;
  Matrix(size_t nelem, size_t nip);
  Matrix(size_t nelem, size_t nip, double kappa, double sig0, double eps0, double m);

  // return shape
  size_t nelem() const;
  size_t nip() const;

  // parameters
  xt::xtensor<double,2> kappa() const;
  xt::xtensor<double,2> sig0() const;
  xt::xtensor<double,2> eps0() const;
  xt::xtensor<double,2> m() const;

  // check that a type has been set everywhere
  void check() const;

  // set parameters
  void set(const xt::xtensor<size_t,2> &I, double kappa, double sig0, double eps0, double m);

  // return unit matrix of unit tensors
  xt::xtensor<double,4> I() const;
  xt::xtensor<double,6> II() const;
  xt::xtensor<double,6> I4() const;
  xt::xtensor<double,6> I4rt() const;
  xt::xtensor<double,6> I4s() const;
  xt::xtensor<double,6> I4d() const;

  // compute stress (no allocation)
  void Sig(const xt::xtensor<double,4> &Eps,
    xt::xtensor<double,4> &Sig) const;

  // compute stress & tangent (no allocation)
  void Tangent(const xt::xtensor<double,4> &Eps,
    xt::xtensor<double,4> &Sig, xt::xtensor<double,6> &Tangent) const;

  // compute stress (return allocated result)
  xt::xtensor<double,4> Sig(const xt::xtensor<double,4> &Eps) const;

  // compute stress & tangent (return allocated result)
  std::tuple<xt::xtensor<double,4>,xt::xtensor<double,6>> Tangent(const xt::xtensor<double,4> &Eps) const;

private:

  xt::xtensor<int   ,2> m_set;
  xt::xtensor<double,2> m_kappa;
  xt::xtensor<double,2> m_sig0;
  xt::xtensor<double,2> m_eps0;
  xt::xtensor<double,2> m_m;
  size_t                m_nelem;
  size_t                m_nip;
  static const size_t   m_ndim=3;
};

// -------------------------------------------------------------------------------------------------

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.hpp"

// -------------------------------------------------------------------------------------------------

#endif

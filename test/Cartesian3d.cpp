
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

#include <GMatNonLinearElastic/Cartesian3d.h>

namespace GM = GMatNonLinearElastic::Cartesian3d;

TEST_CASE("GMatNonLinearElastic::Cartesian3d", "Cartesian3d.h")
{

double kappa = 12.3;
double sig0 = 45.6;
double eps0 = 45.6;
double m = 1.0;

GM::Tensor2 Eps;
GM::Tensor2 Sig;

double epsm = 0.12;

Eps.fill(0.0);
Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm;

SECTION("NonLinearElastic")
{
    GM::NonLinearElastic mat(kappa, sig0, eps0, m);

    Sig = mat.Stress(Eps);

    EQ(Sig(0,0), 3.0 * kappa * epsm);
    EQ(Sig(1,1), 3.0 * kappa * epsm);
    EQ(Sig(2,2), 3.0 * kappa * epsm);
    EQ(Sig(0,1), 0.0);
    EQ(Sig(1,0), 0.0);
    EQ(Sig(0,2), 0.0);
    EQ(Sig(1,2), 0.0);
    EQ(Sig(2,0), 0.0);
    EQ(Sig(2,1), 0.0);
}

SECTION("Matrix")
{
    size_t nelem = 3;
    size_t nip = 2;

    GM::Matrix mat(nelem, nip);

    // all rows elastic
    {
        xt::xtensor<size_t,2> I = xt::ones<size_t>({nelem, nip});
        mat.setNonLinearElastic(I, kappa, sig0, eps0, m);
    }

    xt::xtensor<double,4> eps = xt::empty<double>({nelem, nip, 3ul, 3ul});

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            xt::view(eps, e, q) = Eps;
        }
    }

    auto sig = mat.Stress(eps);

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            EQ(sig(e,q,0,0), 3.0 * kappa * epsm);
            EQ(sig(e,q,1,1), 3.0 * kappa * epsm);
            EQ(sig(e,q,2,2), 3.0 * kappa * epsm);
            EQ(sig(e,q,0,1), 0.0);
            EQ(sig(e,q,1,0), 0.0);
        }
    }

    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(0), xt::keep(2)), 0.0));
    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(1), xt::keep(2)), 0.0));
    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(2), xt::keep(0)), 0.0));
    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::keep(2), xt::keep(1)), 0.0));
}

}

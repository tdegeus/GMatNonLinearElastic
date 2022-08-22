/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATNONLINEARELASTIC_CARTESIAN3D_H
#define GMATNONLINEARELASTIC_CARTESIAN3D_H

#include <GMatElastic/Cartesian3d.h>
#include <GMatTensor/Cartesian3d.h>

#include "config.h"
#include "version.h"

namespace GMatNonLinearElastic {

/**
Implementation in a 3-d Cartesian coordinate frame.
*/
namespace Cartesian3d {

using GMatElastic::Cartesian3d::epseq;
using GMatElastic::Cartesian3d::Epseq;
using GMatElastic::Cartesian3d::sigeq;
using GMatElastic::Cartesian3d::Sigeq;

/**
Array of material points with a linear elasto-plastic constitutive response with linear hardening.

\tparam N Rank of the array.
*/
template <size_t N>
class NonLinearElastic : public GMatElastic::Cartesian3d::Elastic<N> {
private:
    array_type::tensor<double, N> m_kappa; ///< Bulk modulus per item.
    array_type::tensor<double, N> m_sig0; ///< Reference stress per item.
    array_type::tensor<double, N> m_eps0; ///< Reference strain per item.
    array_type::tensor<double, N> m_m; ///< Exponent per item.

    using GMatElastic::Cartesian3d::Elastic<N>::m_Eps;
    using GMatElastic::Cartesian3d::Elastic<N>::m_Sig;
    using GMatElastic::Cartesian3d::Elastic<N>::m_C;

    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;

    using GMatElastic::Cartesian3d::Elastic<N>::energy;
    using GMatElastic::Cartesian3d::Elastic<N>::K;
    using GMatElastic::Cartesian3d::Elastic<N>::G;

public:
    using GMatElastic::Cartesian3d::Elastic<N>::rank;

    NonLinearElastic() = default;

    /**
    Construct system.
    \param kappa Bulk modulus per item.
    \param sig0 Reference stress per item.
    \param eps0 Reference strain per item.
    \param m Exponent per item.
    */
    template <class T>
    NonLinearElastic(const T& kappa, const T& sig0, const T& eps0, const T& m)
    {
        GMATELASTIC_ASSERT(xt::has_shape(kappa, sig0.shape()));
        GMATELASTIC_ASSERT(xt::has_shape(kappa, eps0.shape()));
        GMATELASTIC_ASSERT(xt::has_shape(kappa, m.shape()));

        // allocating parent class
        std::copy(kappa.shape().cbegin(), kappa.shape().cend(), m_shape.begin());
        this->init(m_shape);

        m_kappa = kappa;
        m_sig0 = sig0;
        m_eps0 = eps0;
        m_m = m;

        m_Eps = xt::zeros<double>(m_shape_tensor2);
        m_Sig = xt::empty<double>(m_shape_tensor2);
        m_C = xt::empty<double>(m_shape_tensor4);

        this->refresh(true); // initialise tangent
    }

    /**
    Bulk modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& kappa() const
    {
        return m_kappa;
    }

    /**
    Shear modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& sig0() const
    {
        return m_sig0;
    }

    /**
    Initial yield stress per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& eps0() const
    {
        return m_eps0;
    }

    /**
    Hardening modulus per item.
    \return [shape()].
    */
    const array_type::tensor<double, N>& m() const
    {
        return m_m;
    }

    void refresh(bool compute_tangent = true) override
    {
        namespace GT = GMatTensor::Cartesian3d::pointer;

#pragma omp parallel
        {
            auto II = GMatTensor::Cartesian3d::II();
            auto I4d = GMatTensor::Cartesian3d::I4d();

            double kappa;
            double sig0;
            double eps0;
            double m;

            auto Eps = xt::adapt(m_Eps.data(), {m_ndim, m_ndim});
            auto Sig = xt::adapt(m_Sig.data(), {m_ndim, m_ndim});
            auto C = xt::adapt(m_C.data(), {m_ndim, m_ndim, m_ndim, m_ndim});

            std::array<double, 9> Epsd;

#pragma omp for
            for (size_t i = 0; i < m_size; ++i) {

                kappa = m_kappa.flat(i);
                sig0 = m_sig0.flat(i);
                eps0 = m_eps0.flat(i);
                m = m_m.flat(i);

                Eps.reset_buffer(&m_Eps.flat(i * m_stride_tensor2), m_stride_tensor2);
                Sig.reset_buffer(&m_Sig.flat(i * m_stride_tensor2), m_stride_tensor2);
                C.reset_buffer(&m_C.flat(i * m_stride_tensor4), m_stride_tensor4);

                double epsm = GT::Hydrostatic_deviatoric(Eps.data(), &Epsd[0]);
                double epseq = std::sqrt(2.0 / 3.0 * GT::A2s_ddot_B2s(&Epsd[0], &Epsd[0]));

                if (epseq != 0.0) {
                    double f = 2.0 / 3.0 * sig0 / std::pow(eps0, m) * std::pow(epseq, m - 1.0);
                    Sig.flat(0) = f * Epsd[0] + 3.0 * kappa * epsm;
                    Sig.flat(1) = f * Epsd[1];
                    Sig.flat(2) = f * Epsd[2];
                    Sig.flat(3) = f * Epsd[3];
                    Sig.flat(4) = f * Epsd[4] + 3.0 * kappa * epsm;
                    Sig.flat(5) = f * Epsd[5];
                    Sig.flat(6) = f * Epsd[6];
                    Sig.flat(7) = f * Epsd[7];
                    Sig.flat(8) = f * Epsd[8] + 3.0 * kappa * epsm;
                }
                else {
                    Sig.flat(0) = 3.0 * kappa * epsm;
                    Sig.flat(1) = 0.0;
                    Sig.flat(2) = 0.0;
                    Sig.flat(3) = 0.0;
                    Sig.flat(4) = 3.0 * kappa * epsm;
                    Sig.flat(5) = 0.0;
                    Sig.flat(6) = 0.0;
                    Sig.flat(7) = 0.0;
                    Sig.flat(8) = 3.0 * kappa * epsm;
                }

                if (!compute_tangent) {
                    return;
                }

                if (epseq != 0.0) {
                    GT::A2_dyadic_B2(&Epsd[0], &Epsd[0], C.data());
                    C *= 2.0 / 3.0 * (m - 1.0) * std::pow(epseq, m - 3.0);
                    C += std::pow(epseq, m - 1.0) * I4d;
                    C *= 2.0 / 3.0 * sig0 / std::pow(eps0, m);
                    C += kappa * II;
                }
                else {
                    xt::noalias(C) = kappa * II + I4d;
                }
            }
        }
    }
};

} // namespace Cartesian3d
} // namespace GMatNonLinearElastic

#endif

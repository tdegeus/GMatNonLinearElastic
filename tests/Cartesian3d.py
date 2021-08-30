import unittest
import numpy as np
import GMatNonLinearElastic.Cartesian3d as GMat

class Test_main(unittest.TestCase):

    def test_NonLinearElastic(self):

        kappa = 12.3
        sig0 = 45.6
        eps0 = 45.6
        m = 1.0

        epsm = 0.12

        Eps = np.array(
            [[epsm, 0.0, 0.0],
             [0.0, epsm, 0.0],
             [0.0, 0.0, epsm]])

        Sig = np.array(
            [[3.0 * kappa * epsm, 0.0, 0.0],
             [0.0, 3.0 * kappa * epsm, 0.0],
             [0.0, 0.0, 3.0 * kappa * epsm]])

        mat = GMat.NonLinearElastic(kappa, sig0, eps0, m)
        mat.setStrain(Eps)

        self.assertTrue(np.allclose(mat.Stress(), Sig))

    def test_Array2d(self):

        kappa = 12.3
        sig0 = 45.6
        eps0 = 45.6
        m = 1.0

        epsm = 0.12

        Eps = np.array(
            [[epsm, 0.0, 0.0],
             [0.0, epsm, 0.0],
             [0.0, 0.0, epsm]])

        Sig = np.array(
            [[3.0 * kappa * epsm, 0.0, 0.0],
             [0.0, 3.0 * kappa * epsm, 0.0],
             [0.0, 0.0, 3.0 * kappa * epsm]])

        nelem = 3
        nip = 2
        mat = GMat.Array2d([nelem, nip])
        ndim = 3

        I = np.ones([nelem, nip], dtype='int')
        mat.setNonLinearElastic(I, kappa, sig0, eps0, m)

        eps = np.zeros((nelem, nip, ndim, ndim))
        sig = np.zeros((nelem, nip, ndim, ndim))

        for e in range(nelem):
            for q in range(nip):
                fac = float((e + 1) * nip + (q + 1))
                eps[e, q, :, :] = fac * Eps
                sig[e, q, :, :] = fac * Sig

        mat.setStrain(eps)

        self.assertTrue(np.allclose(mat.Stress(), sig))

if __name__ == '__main__':

    unittest.main()

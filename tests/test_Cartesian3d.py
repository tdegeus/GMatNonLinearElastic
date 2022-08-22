import unittest

import GMatNonLinearElastic.Cartesian3d as GMat
import numpy as np


class Test_main(unittest.TestCase):
    def test_NonLinearElastic(self):

        shape = [3, 4]
        mat = GMat.NonLinearElastic2d(
            kappa=np.random.random(shape),
            sig0=np.random.random(shape),
            eps0=np.random.random(shape),
            m=np.ones(shape),
        )

        epsm = 0.12
        mat.Eps[..., 0, 0] = epsm
        mat.Eps[..., 1, 1] = epsm
        mat.Eps[..., 2, 2] = epsm
        mat.refresh()

        Sig = np.zeros_like(mat.Eps)
        Sig[..., 0, 0] = 3 * mat.kappa * epsm
        Sig[..., 1, 1] = 3 * mat.kappa * epsm
        Sig[..., 2, 2] = 3 * mat.kappa * epsm

        self.assertTrue(np.allclose(mat.Sig, Sig))


if __name__ == "__main__":

    unittest.main()

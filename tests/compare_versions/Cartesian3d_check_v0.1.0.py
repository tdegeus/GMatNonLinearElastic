import h5py
import numpy as np
import GMatNonLinearElastic.Cartesian3d as GMat
import unittest

class Test(unittest.TestCase):

    def test_main(self):

        with h5py.File('Cartesian3d_random.hdf5') as data:

            kappa = data["kappa"][...]
            sig0 = data["sig0"][...]
            eps0 = data["eps0"][...]
            m = data["m"][...]

            mat = GMat.Matrix(kappa.shape[0], kappa.shape[1])

            for i in range(kappa.shape[0]):
                for j in range(kappa.shape[1]):
                    iden = np.zeros(kappa.shape, dtype=bool)
                    iden[i, j] = True
                    mat.setNonLinearElastic(iden, kappa[i, j], sig0[i, j], eps0[i, j], m[i, j])

            for i in range(20):

                                Eps = data[f"/data/{i:d}/Eps"][...]

                self.assertTrue(np.allclose(mat.Stress(Eps), data[f"/data/{i:d}/Stress"][...]))
                self.assertTrue(np.allclose(mat.Tangent(Eps)[1], data[f"/data/{i:d}/Tangent"][...]))

if __name__ == '__main__':

    unittest.main()

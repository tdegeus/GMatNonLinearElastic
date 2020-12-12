import h5py
import numpy as np
import GMatNonLinearElastic.Cartesian3d as GMat

with h5py.File('Cartesian3d_random.hdf5', 'w') as data:

    nelem = 1000
    nip = 4

    shape = np.array([nelem, nip], np.int)

    data['/shape'] = shape

    mat = GMat.Array2d(shape)

    I = np.ones(shape).astype(np.int)
    n = I.size
    kappa = 12.3
    sig0 = 45.6
    eps0 = 45.6
    m = 1.0

    data['/model/I'] = I
    data['/model/kappa'] = kappa
    data['/model/sig0'] = sig0
    data['/model/eps0'] = eps0
    data['/model/m'] = m

    mat.setNonLinearElastic(I, kappa, sig0, eps0, m)

    for i in range(20):

        GradU = 200 * np.random.random([nelem, nip, 3, 3])

        data['/random/{0:d}/GradU'.format(i)] = GradU

        Eps = np.einsum('...ijkl,...lk->...ij', mat.I4s(), GradU)
        mat.setStrain(Eps)

        data['/random/{0:d}/Stress'.format(i)] = mat.Stress()
        data['/random/{0:d}/Tangent'.format(i)] = mat.Tangent()


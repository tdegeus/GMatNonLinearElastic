import GMatNonLinearElastic.Cartesian3d as GMat
import GMatTensor.Cartesian3d as tensor
import h5py
import numpy as np

with h5py.File("Cartesian3d_random.hdf5", "w") as data:

    shape = [1000, 4]

    mat = GMat.NonLinearElastic2d(
        kappa=np.random.random(shape),
        sig0=np.random.random(shape),
        eps0=np.random.random(shape),
        m=0.1 + 0.1 * np.random.random(shape),
    )

    data["kappa"] = mat.kappa
    data["sig0"] = mat.sig0
    data["eps0"] = mat.eps0
    data["m"] = mat.m

    for i in range(20):

        mat.Eps = tensor.Sym(20 * np.random.random(shape + [3, 3]))

        data[f"/data/{i:d}/Eps"] = mat.Eps
        data[f"/data/{i:d}/Stress"] = mat.Sig
        data[f"/data/{i:d}/Tangent"] = mat.C

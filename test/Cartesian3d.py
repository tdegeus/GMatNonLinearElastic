
import GMatNonLinearElastic.Cartesian3d as GMat
import numpy as np

# ==================================================================================================

def EQ(a,b):
  assert np.abs(a-b) < 1.e-12

def ALLEQ(a, b):
  assert np.allclose(a, b)

# ==================================================================================================

# material model
# - parameters
kappa = 12.3
sig0 = 45.6
eps0 = 45.6
m = 1.0
# - model
mat = GMat.NonLinearElastic(kappa, sig0, eps0, m)

# simple shear + volumetric deformation
# - parameters
epsm = 0.12
# - strain
Eps = [[epsm, 0.0 , 0.0 ],
       [0.0 , epsm, 0.0 ],
       [0.0 , 0.0 , epsm]]
# - stress
Sig = mat.Stress(Eps)
# - analytical solution
EQ(Sig[0,0], 3.0 * kappa * epsm)
EQ(Sig[1,1], 3.0 * kappa * epsm)
EQ(Sig[2,2], 3.0 * kappa * epsm)
EQ(Sig[0,1], 0)
EQ(Sig[1,0], 0)
EQ(Sig[0,2], 0)
EQ(Sig[1,2], 0)
EQ(Sig[2,0], 0)
EQ(Sig[2,1], 0)

# ==================================================================================================

# parameters
kappa = 12.3
sig0 = 45.6
eps0 = 45.6
m = 1.0

# allocate matrix
nelem = 2
nip = 2
mat = GMat.Matrix(nelem, nip)

# all rows: elastic
I = np.ones([nelem, nip], dtype='int')
mat.setNonLinearElastic(I, kappa, sig0, eps0, m)

# simple shear + volumetric deformation
# - parameters
epsm = 0.12;
# - strain
Eps = np.zeros((nelem, nip, 3, 3))
Eps[:,:,0,0] = epsm
Eps[:,:,1,1] = epsm
Eps[:,:,2,2] = epsm
# - stress
Sig = mat.Stress(Eps)

# - analytical solution
EQ(Sig[0,0,0,0], 3.0 * kappa * epsm); EQ(Sig[0,1,0,0], 3.0 * kappa * epsm)
EQ(Sig[0,0,1,1], 3.0 * kappa * epsm); EQ(Sig[0,1,1,1], 3.0 * kappa * epsm)
EQ(Sig[0,0,2,2], 3.0 * kappa * epsm); EQ(Sig[0,1,2,2], 3.0 * kappa * epsm)
EQ(Sig[0,0,0,1], 0);                  EQ(Sig[0,1,0,1], 0)
EQ(Sig[0,0,0,1], 0);                  EQ(Sig[0,1,1,0], 0)
EQ(Sig[1,0,0,0], 3.0 * kappa * epsm); EQ(Sig[1,1,0,0], 3.0 * kappa * epsm)
EQ(Sig[1,0,1,1], 3.0 * kappa * epsm); EQ(Sig[1,1,1,1], 3.0 * kappa * epsm)
EQ(Sig[1,0,2,2], 3.0 * kappa * epsm); EQ(Sig[1,1,2,2], 3.0 * kappa * epsm)
EQ(Sig[1,0,0,1], 0);                  EQ(Sig[1,1,0,1], 0)
EQ(Sig[1,0,0,1], 0);                  EQ(Sig[1,1,1,0], 0)
ALLEQ(Sig[:,:,0,2], 0)
ALLEQ(Sig[:,:,1,2], 0)
ALLEQ(Sig[:,:,2,0], 0)
ALLEQ(Sig[:,:,2,1], 0)

# ==================================================================================================

print('All checks passed')

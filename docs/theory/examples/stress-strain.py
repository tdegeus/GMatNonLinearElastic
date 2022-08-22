import GMatNonLinearElastic.Cartesian3d as GMat
import matplotlib.pyplot as plt
import numpy as np

try:
    plt.style.use(["goose", "goose-latex"])
except FileNotFoundError:
    pass

# get response

mat = GMat.NonLinearElastic0d(10.0, 1.0, 1.0, 0.1)

ninc = 301
epseq = np.zeros(ninc)
sigeq = np.zeros(ninc)

for igamma, gamma in enumerate(np.linspace(0.0, 0.1, ninc)):

    mat.Eps = np.array(
        [
            [0.0, gamma, 0.0],
            [gamma, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    epseq[igamma] = GMat.Epseq(mat.Eps)
    sigeq[igamma] = GMat.Sigeq(mat.Sig)

# plot result

fig, ax = plt.subplots()

ax.plot(epseq, sigeq)

ax.set_xlabel(r"$\varepsilon_\mathrm{eq}$")
ax.set_ylabel(r"$\sigma_\mathrm{eq}$")

fig.savefig("stress-strain.pdf")
plt.show()

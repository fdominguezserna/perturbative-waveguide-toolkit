"""Example asymmetric slab calculation."""

import numpy as np
import matplotlib.pyplot as plt

from perturbative_waveguides import find_effective_indices, slab_te_tm_fields_asym

n1 = 1.44
n2 = 2.00
n3 = 1.00
lambda0 = 1.55
d = 0.7
mode = "TE"

neff_grid = np.linspace(min(n1, n2, n3), max(n1, n2, n3), 3500)
neffs = find_effective_indices(d, lambda0, n1, n2, n3, neff_grid, mode)
print("Supported neffs:", neffs)

profile = slab_te_tm_fields_asym(d, lambda0, neffs[0], n1, n2, n3, mode)
print("Core fraction:", profile.core_fraction)

plt.plot(profile.x, abs(profile.field) ** 2)
plt.xlabel("x")
plt.ylabel(r"$|U(x)|^2$")
plt.title("Asymmetric slab field profile")
plt.tight_layout()
plt.show()

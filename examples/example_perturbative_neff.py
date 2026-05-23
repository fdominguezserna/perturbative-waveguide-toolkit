"""Example perturbative effective-index calculation."""

import matplotlib.pyplot as plt

from perturbative_waveguides import effective_index_perturbative



result = effective_index_perturbative(
    n1=1.44,
    n2=2.00,
    n3=1.00,
    lambda0=1.55,
    width=0.8,
    height=0.7,
    mode="TE",
    l=1,
    m=1,
    nps=3500,
)

print(f"Perturbative neff = {result.neff_mode:.8f}")
print(f"Equivalent cladding index = {result.n_equiv_cladding:.8f}")
print("Lateral slab neffs:", result.neffs)

# plt.plot(result.x, abs(result.field) ** 2)
# plt.xlabel("x")
# plt.ylabel(r"$|U(x)|^2$")
# plt.title("Perturbative lateral field profile")
# plt.tight_layout()
# plt.show()

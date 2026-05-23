"""Compare perturbative effective-index calculation against Lumerical data.

Validation Disclaimer

The perturbative model implemented in this repository assumes a silicon nitride (Si₃N₄: 2um x 1um, w x h) waveguide on a silicon dioxide (SiO₂) substrate, using the Sellmeier dispersion coefficients specified below.
The same material dispersion relations were employed in independent simulations performed with Ansys Lumerical MODE Solutions for an equivalent waveguide geometry. The simulation results were exported to:

guia_lum_ajustado_Sellmeier_TM_h1p0_w2p0_100pts.mat

and are included solely as a reference dataset for validating the perturbative calculations against full-wave numerical simulations.

It is intentionally kept as an example/validation script rather than as part of
the main package API.

Run from the project root after installing the package in editable mode:

    pip install -e .
    python examples/compare_lumerical_sellmeier.py

The script loads:
- examples/data/guia_lum_ajustado_Sellmeier_TM_h1p0_w2p0_100pts.mat
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline

from perturbative_waveguides import effective_index_perturbative


C_UM_PER_S = 3.0e14  # speed of light in um/s, matching the MATLAB script


def sellmeier_index(coeffs: np.ndarray | list[float], wavelength_um: np.ndarray | float, model: int) -> np.ndarray:
    """Evaluate the adjusted Sellmeier model.

    Parameters
    ----------
    coeffs:
        Sellmeier coefficients.
    wavelength_um:
        Wavelength in micrometers.
    model:
        ``2`` for a two-resonance model, ``3`` for a three-resonance model
        with the constant offset 11.6 used in the original script.

    Returns
    -------
    ndarray
        Refractive index evaluated at ``wavelength_um``.
    """

    x = np.asarray(coeffs, dtype=float)
    lam = np.asarray(wavelength_um, dtype=float)

    if model == 2:
        n2 = 1.0 + (x[0] * lam**2) / (lam**2 - x[1] ** 2) + (x[2] * lam**2) / (lam**2 - x[3] ** 2)
    elif model == 3:
        n2 = (
            11.6
            + (x[0] * lam**2) / (lam**2 - x[1] ** 2)
            + (x[2] * lam**2) / (lam**2 - x[3] ** 2)
            + (x[4] * lam**2) / (lam**2 - x[5] ** 2)
        )
    else:
        raise ValueError("model must be 2 or 3.")

    return np.sqrt(n2)


def main() -> None:
    data_dir = Path(__file__).resolve().parent / "data"

    lum = loadmat(data_dir / "guia_lum_ajustado_Sellmeier_TM_h1p0_w2p0_100pts.mat")
    f_hz = np.ravel(lum["f"]).astype(float)
    neff_lum = np.real(np.ravel(lum["neff"]))

    

    # Adjusted Sellmeier coefficients for materials dispersions used in numerical simulations.
    x_si = np.array([-8.08431312740552, -3.94071108406619e-05, 8.22519274277156, 0.336004326525049, 8.00161577548198, 279.999907646122])
    x_sio2 = np.array([1.10378795747464, 0.0899411578874705, 0.666874465706605, 8.61718795571091])
    x_si3n4 = np.array([2.97155143840414, 0.135376454869154, 3.10950618314236, 10.5418253140356])

    # Geometry and mode selection.
    wavelengths_um = np.linspace(0.5, 2.0, 100)
    width_um = 2.0
    height_um = 1.0
    l_vertical = 1  # MATLAB-style one-based index
    m_horizontal = 1  # MATLAB-style one-based index
    mode = "TM"
    nps = 8000

    neff_pert = np.zeros_like(wavelengths_um, dtype=float)

    for i, lam in enumerate(wavelengths_um):
        n1 = float(sellmeier_index(x_sio2, lam, 2))
        n2 = float(sellmeier_index(x_si3n4, lam, 2))
        n3 = 1.0

        result = effective_index_perturbative(
            n1=n1,
            n2=n2,
            n3=n3,
            lambda0=lam,
            width=width_um,
            height=height_um,
            mode=mode,
            l=l_vertical,
            m=m_horizontal,
            nps=nps,
        )
        neff_pert[i] = result.neff_mode

    lambda_lum_um = C_UM_PER_S / f_hz

    fig1, ax1 = plt.subplots()
    ax1.plot(lambda_lum_um, neff_lum, "ro", label="Lumerical")
    ax1.plot(wavelengths_um, neff_pert, "-k", label="Perturbative")
    ax1.set_xlabel(r"$\lambda$ ($\mu$m)")
    ax1.set_ylabel(r"$n_\mathrm{eff}$")
    ax1.legend()
    ax1.tick_params(direction="out")

    # Compare propagation constants k = 2*pi*n_eff/lambda.
    k_lum = (2.0 * np.pi / lambda_lum_um) * neff_lum
    sort_idx = np.argsort(lambda_lum_um)
    spline_neff = CubicSpline(lambda_lum_um[sort_idx], neff_lum[sort_idx])
    neff_lum_interp = spline_neff(wavelengths_um)
    k_lum_interp = (2.0 * np.pi / wavelengths_um) * neff_lum_interp
    k_pert = (2.0 * np.pi / wavelengths_um) * neff_pert


    with np.errstate(all="ignore"):
        coeff_lum = np.polyfit(wavelengths_um, k_lum_interp, 35)

    fig2, ax2 = plt.subplots()
    ax2.plot(lambda_lum_um, k_lum, "*r", label="Lumerical")
    ax2.plot(wavelengths_um, k_pert, "-k", label="Perturbative")
    ax2.set_xlabel(r"$\lambda$ ($\mu$m)")
    ax2.set_ylabel(r"$k$ ($1/\mu$m)")
    ax2.legend()
    ax2.tick_params(direction="out")

    # A compact numerical check useful when running from the terminal.
    rms_neff = np.sqrt(np.mean((neff_lum_interp - neff_pert) ** 2))
    max_abs_neff = np.max(np.abs(neff_lum_interp - neff_pert))
    print(f"RMS difference in neff: {rms_neff:.6e}")
    print(f"Max absolute difference in neff: {max_abs_neff:.6e}")
    #print(f"Polynomial coefficients for interpolated Lumerical k(lambda), degree 35:")
    #print(coeff_lum)

    plt.show()


if __name__ == "__main__":
    main()

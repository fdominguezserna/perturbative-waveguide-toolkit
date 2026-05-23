"""Slab-waveguide mode equations and field profiles.

The implementation follows the notation of the original scripts:

- ``n1``: lower/left cladding index
- ``n2``: core index
- ``n3``: upper/right cladding index
- ``d``: slab thickness/width
- ``lambda0``: vacuum wavelength, in the same length units as ``d``
- ``mode``: ``"TE"`` or ``"TM"``

All lengths must use a consistent unit system, for example micrometers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import root_scalar

Mode = Literal["TE", "TM"]


@dataclass(frozen=True)
class SlabFieldProfile:
    """Field profile returned by :func:`slab_te_tm_fields_asym`.

    Attributes
    ----------
    field:
        Complex field profile sampled on ``x``.
    x:
        Spatial grid.
    core_fraction:
        Fraction of ``|field|^2`` contained in the slab core ``0 < x < d``.
    """

    field: np.ndarray
    x: np.ndarray
    core_fraction: float


def _validate_mode(mode: str) -> Mode:
    mode_upper = mode.upper()
    if mode_upper not in {"TE", "TM"}:
        raise ValueError("mode must be 'TE' or 'TM'.")
    return mode_upper  # type: ignore[return-value]


def _as_real_if_valid(value: np.ndarray | complex | float) -> np.ndarray:
    """Return real part where imaginary residue is numerical noise.
    """

    arr = np.asarray(value, dtype=np.complex128)
    out = np.empty(arr.shape, dtype=float)
    valid = np.isclose(arr.imag, 0.0, atol=1e-10, rtol=1e-10)
    out[valid] = arr.real[valid]
    out[~valid] = 1.0e6
    return out


def characteristic_difference(
    d: float,
    lambda0: float,
    n1: float,
    n2: float,
    n3: float,
    neff: float | np.ndarray,
    mode: str = "TE",
) -> float | np.ndarray:
    """Evaluate the asymmetric-slab characteristic equation residual.

    This is the Python equivalent of ``neffA.m``. Its roots give the supported
    effective indices of an asymmetric slab waveguide.

    Parameters
    ----------
    d, lambda0:
        Slab thickness and vacuum wavelength, using consistent length units.
    n1, n2, n3:
        Refractive indices of the two claddings and the core.
    neff:
        Trial effective index or array of trial effective indices.
    mode:
        ``"TE"`` or ``"TM"``.

    Returns
    -------
    float or ndarray
        Real residual ``A - B``. Non-real values are replaced by ``1e6``
    """

    mode = _validate_mode(mode)
    neff_arr = np.asarray(neff, dtype=np.complex128)

    k0 = 2.0 * np.pi / lambda0
    beta = k0 * neff_arr
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        kx = np.sqrt((k0 * n2) ** 2 - beta**2)

        a_term = np.tan(d * k0 * np.sqrt(n2**2 - neff_arr**2))
        gamma1 = np.sqrt(beta**2 - (k0 * n1) ** 2)
        gamma3 = np.sqrt(beta**2 - (k0 * n3) ** 2)

        if mode == "TE":
            b_term = ((gamma1 / kx) + (gamma3 / kx)) / (1.0 - gamma1 * gamma3 / kx**2)
        else:
            b_term = (
                (n2 / n1) ** 2 * (gamma1 / kx)
                + (n2 / n3) ** 2 * (gamma3 / kx)
            ) / (1.0 - gamma1 * gamma3 * (n2**4 / (n1**2 * n3**2)) / kx**2)

        residual = _as_real_if_valid(a_term - b_term)
    if np.ndim(neff) == 0:
        return float(np.asarray(residual))
    return residual


def find_effective_indices(
    d: float,
    lambda0: float,
    n1: float,
    n2: float,
    n3: float,
    neff_grid: np.ndarray,
    mode: str = "TE",
    *,
    residual_tol: float = 1.0e-10,
    xtol: float = 1.0e-14,
    rtol: float = 1.0e-14,
) -> np.ndarray:
    """Find slab-waveguide effective indices by bracketing sign changes.

    Parameters
    ----------
    d, lambda0, n1, n2, n3, mode:
        Waveguide parameters. See :func:`characteristic_difference`.
    neff_grid:
        One-dimensional array of trial effective indices. A dense grid gives
        more reliable root detection near tangent discontinuities.
    residual_tol:
        Accepted absolute residual after root refinement.
    xtol, rtol:
        Absolute and relative tolerances passed to SciPy's root finder.

    Returns
    -------
    ndarray
        Effective indices sorted in descending order.
    """

    neff_grid = np.asarray(neff_grid, dtype=float)
    if neff_grid.ndim != 1 or neff_grid.size < 2:
        raise ValueError("neff_grid must be a one-dimensional array with at least two points.")

    residuals = np.asarray(
        characteristic_difference(d, lambda0, n1, n2, n3, neff_grid, mode),
        dtype=float,
    )

    finite = np.isfinite(residuals)
    candidates: list[float] = []

    for i in range(neff_grid.size - 1):
        if not (finite[i] and finite[i + 1]):
            continue

        f_left = residuals[i]
        f_right = residuals[i + 1]
        if f_left == 0.0:
            candidates.append(neff_grid[i])
            continue
        if f_left * f_right > 0.0:
            continue

        a = neff_grid[i]
        b = neff_grid[i + 1]
        if a == b:
            continue

        try:
            result = root_scalar(
                lambda x: characteristic_difference(d, lambda0, n1, n2, n3, x, mode),
                bracket=(a, b),
                method="brentq",
                xtol=xtol,
                rtol=rtol,
            )
        except ValueError:
            continue

        if result.converged:
            root = float(result.root)
            f_root = float(characteristic_difference(d, lambda0, n1, n2, n3, root, mode))
            if abs(f_root) < residual_tol:
                candidates.append(root)

    if not candidates:
        return np.array([], dtype=float)

    # Remove duplicates caused by grid points falling very close to the same root.
    candidates_sorted = np.sort(np.array(candidates, dtype=float))[::-1]
    unique = []
    for value in candidates_sorted:
        if not unique or not np.isclose(value, unique[-1], rtol=1e-10, atol=1e-12):
            unique.append(value)
    return np.array(unique, dtype=float)


def slab_te_tm_fields_asym(
    d: float,
    lambda0: float,
    neff: float,
    n1: float,
    n2: float,
    n3: float,
    mode: str = "TE",
    *,
    npts: int = 500,
) -> SlabFieldProfile:
    """Compute TE/TM field profile for an asymmetric slab waveguide.
    """

    mode = _validate_mode(mode)
    if npts < 2:
        raise ValueError("npts must be at least 2.")

    k0 = 2.0 * np.pi / lambda0
    beta = k0 * neff
    kx = np.sqrt((k0 * n2) ** 2 - beta**2 + 0j)
    gamma1 = np.sqrt(beta**2 - (k0 * n1) ** 2 + 0j)
    gamma3 = np.sqrt(beta**2 - (k0 * n3) ** 2 + 0j)

    x = np.linspace(-1.5 * d, 2.5 * d, npts)

    b_amp = 1.0
    d_amp = 1.0
    if mode == "TE":
        a_amp = gamma1 / kx
    else:
        a_amp = (n2 / n1) ** 2 * gamma1 / kx

    c_amp = np.exp(gamma3 * d) * (a_amp * np.sin(kx * d) + b_amp * np.cos(kx * d))

    core_mask = (x < d) & (x > 0.0)
    lower_mask = x <= 0.0
    upper_mask = x >= d

    field = (
        (a_amp * np.sin(kx * x) + b_amp * np.cos(kx * x)) * core_mask
        + d_amp * np.exp(gamma1 * x) * lower_mask
        + c_amp * np.exp(-gamma3 * x) * upper_mask
    )

    dx = x[1] - x[0]
    core_power = float(np.sum(np.abs(field * core_mask) ** 2) * dx)
    total_power = float(np.sum(np.abs(field) ** 2) * dx)
    core_fraction = core_power / total_power if total_power != 0.0 else np.nan

    return SlabFieldProfile(field=field, x=x, core_fraction=core_fraction)


# MATLAB-style aliases for users who want a direct mapping from the old scripts.
neffA = characteristic_difference
neffO = find_effective_indices
SlabTEMFields_Asym = slab_te_tm_fields_asym

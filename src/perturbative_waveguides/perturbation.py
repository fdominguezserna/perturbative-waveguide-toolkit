"""Perturbative effective-index workflow for rectangular waveguides.

The method first solves a vertical asymmetric slab, estimates a perturbative cladding/core correction, and
then solves the lateral slab using the corrected equivalent index.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .slab import find_effective_indices, slab_te_tm_fields_asym


@dataclass(frozen=True)
class PerturbativeResult:
    """Result of :func:`effective_index_perturbative`.

    Attributes
    ----------
    neff_mode:
        Selected effective index after the lateral slab calculation.
    neffs:
        All effective indices found in the lateral slab step, sorted descending.
    field:
        Lateral field profile of the selected mode.
    x:
        Lateral spatial grid.
    n_equiv_cladding:
        Perturbatively corrected equivalent cladding index, named ``n_1`` in
        the original MATLAB script.
    vertical_neffs:
        Effective indices found in the initial vertical slab step.
    vertical_core_fraction:
        Fraction of vertical field intensity inside the vertical slab core.
    """

    neff_mode: float
    neffs: np.ndarray
    field: np.ndarray
    x: np.ndarray
    n_equiv_cladding: float
    vertical_neffs: np.ndarray
    vertical_core_fraction: float


def effective_index_perturbative(
    n1: float,
    n2: float,
    n3: float,
    lambda0: float,
    width: float,
    height: float,
    mode: str,
    l: int,
    m: int,
    nps: int,
) -> PerturbativeResult:
    """Compute a perturbative effective index for a rectangular waveguide.

    This is the Python equivalent of ``neffPerturb.m``.

    Parameters
    ----------
    n1, n2, n3:
        Refractive indices used in the initial vertical slab calculation.
    lambda0:
        Vacuum wavelength, in the same length units as ``width`` and ``height``.
    width, height:
        Rectangular waveguide dimensions.
    mode:
        ``"TE"`` or ``"TM"``.
    l:
        One-based vertical slab mode number, matching the MATLAB script.
    m:
        One-based lateral slab mode number, matching the MATLAB script.
    nps:
        Number of points in the effective-index search grid.

    Returns
    -------
    PerturbativeResult
        Dataclass containing the selected mode, all lateral modes, lateral
        field profile, and the corrected equivalent cladding index.

    Notes
    -----
    ``l`` and ``m`` are intentionally one-based. For Python-native zero-based mode selection, subtract
    one before calling or wrap this function in your own interface.
    """

    if nps < 2:
        raise ValueError("nps must be at least 2.")
    if l < 1 or m < 1:
        raise ValueError("l and m are one-based mode numbers and must be >= 1.")

    nrange = np.linspace(min(n1, n2, n3), max(n1, n2, n3), nps)
    vertical_neffs = find_effective_indices(height, lambda0, n1, n2, n3, nrange, mode)

    if l > len(vertical_neffs):
        raise ValueError(f"Vertical slab does not support l={l}; found {len(vertical_neffs)} modes.")

    vertical_neff = float(vertical_neffs[l - 1])
    vertical_profile = slab_te_tm_fields_asym(height, lambda0, vertical_neff, n1, n2, n3, mode)
    frac = vertical_profile.core_fraction

    na = n2
    delta_na = -(n2 - n3)
    n_barra = vertical_neff
    n_av = (2.0 * na + delta_na) / 2.0
    delta_n = (n_av * delta_na / n_barra) * frac
    n_equiv_cladding = n_barra + delta_n

    lateral_n1 = n_equiv_cladding
    lateral_n2 = vertical_neff
    lateral_n3 = n_equiv_cladding

    lateral_range = np.linspace(
        min(lateral_n1, lateral_n2, lateral_n3),
        max(lateral_n1, lateral_n2, lateral_n3),
        nps,
    )
    lateral_neffs = find_effective_indices(
        width, lambda0, lateral_n1, lateral_n2, lateral_n3, lateral_range, mode
    )

    if m > len(lateral_neffs):
        raise ValueError(
            f"Lateral slab does not support m={m}; vertical modes={len(vertical_neffs)}, "
            f"lateral modes={len(lateral_neffs)}."
        )

    neff_mode = float(lateral_neffs[m - 1])
    lateral_profile = slab_te_tm_fields_asym(
        width, lambda0, neff_mode, lateral_n1, lateral_n2, lateral_n3, mode
    )

    return PerturbativeResult(
        neff_mode=neff_mode,
        neffs=lateral_neffs,
        field=lateral_profile.field,
        x=lateral_profile.x,
        n_equiv_cladding=float(np.real_if_close(n_equiv_cladding)),
        vertical_neffs=vertical_neffs,
        vertical_core_fraction=float(frac),
    )


# MATLAB-style alias.
neffPerturb = effective_index_perturbative

import numpy as np

from perturbative_waveguides import find_effective_indices, effective_index_perturbative


def test_slab_modes_are_found():
    neff_grid = np.linspace(1.0, 2.0, 1500)
    neffs = find_effective_indices(0.7, 1.55, 1.44, 2.0, 1.0, neff_grid, "TE")
    assert len(neffs) >= 1
    assert np.all(np.diff(neffs) <= 0)


def test_perturbative_result_is_finite():
    result = effective_index_perturbative(1.44, 2.0, 1.0, 1.55, 0.8, 0.7, "TE", 1, 1, 1500)
    assert np.isfinite(result.neff_mode)
    assert len(result.neffs) >= 1

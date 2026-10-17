"""Test the cycle sub-model machinery and its SPECTRA dependency."""

import pytest

from thermoptcobrapy import ThermOptError
from thermoptcobrapy._tic import consistent_reaction_ids

from .conftest import milp


@milp
def test_missing_spectra_cobra_is_a_clear_error(monkeypatch, cycle_model):
    """Test that a missing dependency fails loudly, not by degrading silently.

    SPECTRA is required, not optional, so losing it should look nothing like
    "no cycles found" -- it should say exactly what is missing and how to
    get it.

    """
    monkeypatch.setitem(__import__("sys").modules, "spectra_cobra", None)

    with pytest.raises(ThermOptError, match="spectra_cobra is not installed"):
        consistent_reaction_ids(cycle_model, tol=1e-4)

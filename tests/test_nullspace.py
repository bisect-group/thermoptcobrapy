"""Test Fast-SNP."""

import numpy as np
from cobra.util.array import create_stoichiometric_matrix

from thermoptcobrapy import cycle_basis, fast_snp, subset_model
from thermoptcobrapy._nullspace import nullity
from thermoptcobrapy._tic import cycle_submodel


def test_basis_vectors_lie_in_the_null_space(cycle_model):
    """Test that every vector returned satisfies S v = 0."""
    ids, basis = cycle_basis(cycle_model, tol=1e-4, seed=0)
    closed = subset_model(cycle_submodel(cycle_model), set(ids))
    stoichiometry = np.asarray(create_stoichiometric_matrix(closed, array_type="dense"))

    assert basis.shape[1] > 0
    assert np.abs(stoichiometry @ basis).max() < 1e-9


def test_basis_rows_follow_the_cycle_reactions(cycle_model):
    """Test that the basis is indexed by the reactions it is reported for."""
    ids, basis = cycle_basis(cycle_model, tol=1e-4, seed=0)

    assert ids == ["R1", "R2", "R3"]
    assert basis.shape[0] == len(ids)
    # The loop runs all three reactions together, so no component is zero.
    assert np.all(np.abs(basis[:, 0]) > 0)


def test_an_acyclic_model_has_an_empty_basis(acyclic_model):
    """Test that a model with no cycle yields nothing to constrain."""
    ids, basis = cycle_basis(acyclic_model, tol=1e-4, seed=0)

    assert ids == []
    assert basis.size == 0


def test_the_basis_is_reproducible_given_a_seed(double_cycle_model):
    """Test that the random weights are controlled by the seed."""
    _, first = cycle_basis(double_cycle_model, tol=1e-4, seed=3)
    _, again = cycle_basis(double_cycle_model, tol=1e-4, seed=3)

    assert np.allclose(first, again)


def test_nullity_counts_the_free_directions():
    """Test the nullity helper against a matrix whose rank is known."""
    # Two rows, three columns, rank 2, so one free direction.
    matrix = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]])
    assert nullity(matrix) == 1


def test_fast_snp_stops_when_the_bounds_allow_nothing(acyclic_model):
    """Test that the search ends rather than returning spurious vectors.

    With the exchanges shut an acyclic network can carry no flux at all, so
    Fast-SNP must come back empty even though ``S`` has a null space.

    """
    closed = cycle_submodel(acyclic_model)
    basis = fast_snp(closed, v_tol=1e-6, seed=0)

    assert basis.shape[1] == 0

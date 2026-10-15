"""Test thermodynamically feasible context-specific extraction."""

import pytest

from thermoptcobrapy import BLOCKED, thermopt_cc, thermoptics
from thermoptcobrapy.exceptions import ThermOptInfeasibleCoreError

from .conftest import milp


@milp
def test_the_core_survives_extraction(cycle_model):
    """Test that the reactions asked for are in the result."""
    extracted, blocked_core = thermoptics(cycle_model, ["R4"], tol=1e-4, seed=0)

    assert blocked_core == []
    assert "R4" in {rxn.id for rxn in extracted.reactions}


@milp
def test_the_extracted_model_carries_no_cycle(cycle_model):
    """Test that the result is free of thermodynamically blocked reactions.

    This is what the method is for: a context-specific model built on a
    network that still contains cycles inherits them.

    """
    extracted, _ = thermoptics(cycle_model, ["R4"], tol=1e-4, seed=0)
    verdicts, _ = thermopt_cc(extracted, tol=1e-4, seed=1)

    assert BLOCKED not in set(verdicts.values())
    # The loop cannot support anything, so none of it should be kept.
    assert {"R1", "R2", "R3"} & {rxn.id for rxn in extracted.reactions} == set()


@milp
def test_a_blocked_core_reaction_is_reported_and_dropped(cycle_model):
    """Test that asking for a cycle reaction warns rather than failing."""
    with pytest.warns(UserWarning, match="thermodynamically blocked"):
        extracted, blocked_core = thermoptics(
            cycle_model, ["R1", "R4"], tol=1e-4, seed=0
        )

    assert blocked_core == ["R1"]
    assert "R1" not in {rxn.id for rxn in extracted.reactions}
    assert "R4" in {rxn.id for rxn in extracted.reactions}


@milp
def test_a_blocked_core_without_the_check_is_an_error(cycle_model):
    """Test that skipping the consistency check surfaces the problem.

    With ``do_tocc`` off the blocked reactions are not removed first, so an
    impossible core has to fail rather than quietly return something wrong.

    """
    with pytest.raises(ThermOptInfeasibleCoreError):
        thermoptics(cycle_model, ["R1"], tol=1e-4, do_tocc=False, seed=0)


@milp
def test_an_unknown_core_reaction_is_rejected(cycle_model):
    """Test that a typo in the core list is caught rather than ignored."""
    with pytest.raises(KeyError, match="not in the model"):
        thermoptics(cycle_model, ["NOPE"], tol=1e-4, seed=0)


@milp
def test_the_input_model_is_not_modified(cycle_model):
    """Test that extraction leaves the reference model alone."""
    before = {rxn.id: rxn.bounds for rxn in cycle_model.reactions}
    thermoptics(cycle_model, ["R4"], tol=1e-4, seed=0)

    assert {rxn.id: rxn.bounds for rxn in cycle_model.reactions} == before


@milp
def test_annotations_survive_extraction(cycle_model):
    """Test that the result keeps what the source carried."""
    cycle_model.reactions.R4.annotation = {"kegg.reaction": "R00001"}
    cycle_model.reactions.R4.name = "the useful one"
    extracted, _ = thermoptics(cycle_model, ["R4"], tol=1e-4, seed=0)

    kept = extracted.reactions.get_by_id("R4")
    assert kept.annotation == {"kegg.reaction": "R00001"}
    assert kept.name == "the useful one"

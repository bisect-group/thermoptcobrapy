"""Test the thermodynamically constrained consistency check."""

import pytest
from cobra.flux_analysis import flux_variability_analysis

import thermoptcobrapy.consistency as consistency_module
from thermoptcobrapy import (
    BLOCKED,
    FORWARD,
    REVERSIBLE,
    UNBLOCKED,
    ThermOptTimeLimitError,
    thermodynamically_blocked_ids,
    thermopt_cc,
)

from .conftest import milp


@milp
def test_a_cycle_is_blocked_though_fva_permits_it(cycle_model):
    """Test the whole point of the method.

    ``R1``, ``R2`` and ``R3`` can only carry flux as a closed loop. The mass
    balance allows it, so FVA reports the full range; thermodynamics does
    not, so the check must block all three while leaving the rest alone.

    """
    ranges = flux_variability_analysis(
        cycle_model, fraction_of_optimum=0.0, processes=1
    )
    assert ranges.loc["R1", "maximum"] > 1.0  # FVA sees nothing wrong

    verdicts, working = thermopt_cc(cycle_model, tol=1e-4, seed=0)

    assert verdicts["R1"] == BLOCKED
    assert verdicts["R2"] == BLOCKED
    assert verdicts["R3"] == BLOCKED
    assert verdicts["R4"] == UNBLOCKED
    assert verdicts["EX_A"] == UNBLOCKED


@milp
def test_directions_agree_with_the_plain_verdict(cycle_model):
    """Test that asking for directions does not change what is blocked."""
    plain, _ = thermopt_cc(cycle_model, tol=1e-4, seed=0)
    detailed, _ = thermopt_cc(cycle_model, tol=1e-4, get_direction=True, seed=0)

    for rxn_id, verdict in plain.items():
        if verdict == BLOCKED:
            assert detailed[rxn_id] == BLOCKED
        else:
            assert detailed[rxn_id] != BLOCKED


@milp
def test_a_reversible_usable_reaction_is_reported_reversible(cycle_model):
    """Test that both directions are found when both are feasible."""
    verdicts, _ = thermopt_cc(cycle_model, tol=1e-4, get_direction=True, seed=0)

    assert verdicts["R4"] == REVERSIBLE


@milp
def test_an_irreversible_reaction_is_reported_forward(cycle_model):
    """Test that a one-way reaction is not called reversible."""
    cycle_model.reactions.R4.bounds = (0.0, 1000.0)
    verdicts, _ = thermopt_cc(cycle_model, tol=1e-4, get_direction=True, seed=0)

    assert verdicts["R4"] == FORWARD


@milp
def test_an_acyclic_model_blocks_nothing(acyclic_model):
    """Test that a network without cycles comes back intact."""
    verdicts, _ = thermopt_cc(acyclic_model, tol=1e-4, seed=0)

    assert set(verdicts.values()) == {UNBLOCKED}


@milp
def test_results_refer_to_the_oriented_model(reverse_only_model):
    """Test that the returned model is the one the verdicts are keyed to.

    The reverse-only reaction is flipped, so the model handed back differs
    from the input and is what the caller has to use.

    """
    verdicts, working = thermopt_cc(reverse_only_model, tol=1e-4, seed=0)

    assert set(verdicts) == {rxn.id for rxn in working.reactions}
    assert working.reactions.R1.bounds == (0.0, 1000.0)
    assert reverse_only_model.reactions.R1.bounds == (-1000.0, 0.0)


@milp
def test_the_input_model_is_not_modified(cycle_model):
    """Test that nothing is left behind on the model handed in."""
    before = (
        len(cycle_model.variables),
        len(cycle_model.constraints),
        {rxn.id: rxn.bounds for rxn in cycle_model.reactions},
    )
    thermopt_cc(cycle_model, tol=1e-4, seed=0)

    assert before == (
        len(cycle_model.variables),
        len(cycle_model.constraints),
        {rxn.id: rxn.bounds for rxn in cycle_model.reactions},
    )


@milp
def test_hitting_the_time_limit_ceiling_is_an_error(monkeypatch, cycle_model):
    """Test that an unproven verdict at the ceiling raises rather than warns.

    Reporting the unsettled reactions as blocked would dress a timeout up as
    an answer, so this has to fail loudly instead.

    """

    class _TimedOut:
        status = "time_limit"

        def active_ids(self):
            return []

    monkeypatch.setattr(
        consistency_module, "solve_indicator_problem", lambda *a, **k: _TimedOut()
    )

    with pytest.raises(ThermOptTimeLimitError):
        thermopt_cc(cycle_model, tol=1e-4, seed=0, time_limit=1.0, max_time_limit=1.0)


@milp
def test_blocked_ids_lists_the_blocked_reactions(cycle_model):
    """Test the convenience wrapper."""
    assert thermodynamically_blocked_ids(cycle_model, tol=1e-4, seed=0) == [
        "R1",
        "R2",
        "R3",
    ]

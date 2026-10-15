"""Test the cycle enumerators."""

import numpy as np
import pytest
from cobra.util.array import create_stoichiometric_matrix

from thermoptcobrapy import reaction_order, thermopt_enum_lp, thermopt_enum_milp
from thermoptcobrapy.enumeration import _close_exchanges

from .conftest import milp


def _as_sets(cycles, directions):
    """Return the cycles as (frozenset, sign-tuple) pairs for comparison."""
    out = set()
    for cycle, direction in zip(cycles, directions):
        smallest = min(abs(v) for v in direction)
        out.add(
            tuple(
                sorted(
                    (rxn, round(value / smallest, 6))
                    for rxn, value in zip(cycle, direction)
                )
            )
        )
    return out


@pytest.mark.parametrize("enumerate_cycles", [thermopt_enum_lp, thermopt_enum_milp])
@milp
def test_a_single_cycle_is_found_both_ways(cycle_model, enumerate_cycles):
    """Test that the loop is found running forwards and backwards."""
    cycles, directions, cycle_rxns, _ = enumerate_cycles(cycle_model)

    assert cycle_rxns == ["R1", "R2", "R3"]
    assert _as_sets(cycles, directions) == {
        (("R1", 1.0), ("R2", 1.0), ("R3", 1.0)),
        (("R1", -1.0), ("R2", -1.0), ("R3", -1.0)),
    }


@pytest.mark.parametrize("enumerate_cycles", [thermopt_enum_lp, thermopt_enum_milp])
@milp
def test_overlapping_cycles_are_all_found(double_cycle_model, enumerate_cycles):
    """Test that every cycle of the network is reported, not just a basis.

    With ``R5`` added there are three cycles, each in two directions: the two
    short loops and the composite that runs ``R5`` against the others.

    """
    cycles, directions, _, _ = enumerate_cycles(double_cycle_model)

    assert _as_sets(cycles, directions) == {
        (("R1", 1.0), ("R5", 1.0)),
        (("R1", -1.0), ("R5", -1.0)),
        (("R1", 1.0), ("R2", 1.0), ("R3", 1.0)),
        (("R1", -1.0), ("R2", -1.0), ("R3", -1.0)),
        (("R2", 1.0), ("R3", 1.0), ("R5", -1.0)),
        (("R2", -1.0), ("R3", -1.0), ("R5", 1.0)),
    }


@pytest.mark.parametrize("enumerate_cycles", [thermopt_enum_lp, thermopt_enum_milp])
@milp
def test_an_acyclic_model_yields_nothing(acyclic_model, enumerate_cycles):
    """Test that a network without cycles reports none."""
    cycles, directions, cycle_rxns, _ = enumerate_cycles(acyclic_model)

    assert cycles == [] and directions == [] and cycle_rxns == []


@milp
def test_reported_cycles_are_steady_state_modes(double_cycle_model):
    """Test that each cycle balances every metabolite it touches."""
    cycles, directions, _, working = thermopt_enum_lp(double_cycle_model)
    closed = _close_exchanges(working)
    stoichiometry = np.asarray(create_stoichiometric_matrix(closed, array_type="dense"))
    position = {rxn.id: i for i, rxn in enumerate(closed.reactions)}

    assert cycles
    for cycle, direction in zip(cycles, directions):
        vector = np.zeros(len(closed.reactions))
        for rxn_id, value in zip(cycle, direction):
            vector[position[rxn_id]] = value
        assert np.abs(stoichiometry @ vector).max() < 1e-6


@milp
def test_the_returned_model_keeps_its_exchanges_open(cycle_model):
    """Test that the reported model is the oriented one, not the closed one.

    The cycles are searched for with the exchanges shut, but the model handed
    back is the one the caller works with, which still has them.

    """
    _, _, _, working = thermopt_enum_lp(cycle_model)

    assert working.reactions.EX_A.bounds != (0.0, 0.0)


def test_reaction_order_puts_small_components_first():
    """Test the ordering used to pick which reaction to search from."""
    from cobra import Metabolite, Model, Reaction

    model = Model("two components")
    reactions = []
    # One component of three reactions, one of a single reaction.
    for rxn_id, stoichiometry in (
        ("big1", {"A": -1, "B": 1}),
        ("big2", {"B": -1, "C": 1}),
        ("big3", {"C": -1, "A": 1}),
        ("small1", {"X": -1, "Y": 1}),
    ):
        reaction = Reaction(rxn_id, lower_bound=-1.0, upper_bound=1.0)
        reaction.add_metabolites({Metabolite(k): v for k, v in stoichiometry.items()})
        reactions.append(reaction)
    model.add_reactions(reactions)

    assert reaction_order(model)[0] == "small1"

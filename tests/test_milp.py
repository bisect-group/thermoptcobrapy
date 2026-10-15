"""Test the two formulations of the indicator program.

The constraints can be expressed either with big-M coefficients, which every
solver understands, or as native indicator constraints, which Gurobi enforces
by branching. They must agree, and the native one must not leave anything
behind on the solver.
"""

import numpy as np
import pytest

from thermoptcobrapy import cycle_basis, prepare
from thermoptcobrapy._milp import (
    SAFETY_FACTOR,
    _big_m_resolves,
    check_tolerances,
    solve_indicator_problem,
    supports_indicators,
)
from thermoptcobrapy.exceptions import ThermOptError

from .conftest import milp


@milp
def test_both_formulations_agree(double_cycle_model):
    """Test that big-M and native indicators reach the same answer."""
    working, _ = prepare(double_cycle_model)
    tic_ids, basis = cycle_basis(working, tol=1e-4, seed=0)
    core = {rxn.id for rxn in working.reactions}

    answers = []
    for native in (True, False):
        problem = solve_indicator_problem(
            working,
            tol=1e-4,
            tic_ids=tic_ids,
            basis=basis,
            objective_ids=core,
            maximise=True,
            minimum_for=core,
            time_limit=60,
            use_indicators=native,
        )
        answers.append(sorted(problem.active_ids()))

    assert answers[0] == answers[1]
    # The loops cannot carry flux, so only the rest should be live.
    assert set(answers[0]) == {"R4", "EX_A", "EX_D"}


@milp
def test_native_indicators_are_removed_again(cycle_model):
    """Test that nothing is left on the solver after a native solve.

    The indicator constraints are attached outside optlang's bookkeeping, so
    the context manager cannot clean them up and the code has to.

    """
    working, _ = prepare(cycle_model)
    tic_ids, basis = cycle_basis(working, tol=1e-4, seed=0)
    before = working.solver.problem.NumGenConstrs

    solve_indicator_problem(
        working,
        tol=1e-4,
        tic_ids=tic_ids,
        basis=basis,
        objective_ids={rxn.id for rxn in working.reactions},
        maximise=True,
        time_limit=60,
        use_indicators=True,
    )
    working.solver.problem.update()

    assert working.solver.problem.NumGenConstrs == before


@milp
def test_native_indicators_lift_the_threshold_floor(cycle_model):
    """Test that the native path accepts a threshold big-M cannot.

    With bounds of 1000 and an integrality tolerance of 1e-9, big-M leaves
    1e-6 of slack on every indicator, so a threshold near it means nothing.
    Indicator constraints carry no such coefficient.

    """
    working, _ = prepare(cycle_model)
    tiny = 1e-7

    with pytest.raises(ThermOptError, match="too small for this model"):
        check_tolerances(working, tiny, 1e-9, 1e-9, big_m=True)

    # No complaint for the native path, and it still blocks the cycle.
    check_tolerances(working, tiny, 1e-9, 1e-9, big_m=False)


def test_the_floor_is_computed_from_the_largest_bound(cycle_model):
    """Test that the rejection threshold follows the model's own bounds."""
    for reaction in cycle_model.reactions:
        reaction.bounds = (-10.0, 10.0)

    # 10 * 1e-9 of slack, so the floor is SAFETY_FACTOR times that.
    floor = SAFETY_FACTOR * 10.0 * 1e-9
    check_tolerances(cycle_model, floor, 1e-9, 1e-12, big_m=True)
    with pytest.raises(ThermOptError, match="too small"):
        check_tolerances(cycle_model, floor / 10.0, 1e-9, 1e-12, big_m=True)


def test_a_threshold_under_the_feasibility_tolerance_is_rejected(cycle_model):
    """Test the second floor, which applies to both formulations."""
    with pytest.raises(ThermOptError, match="feasibility tolerance"):
        check_tolerances(cycle_model, 1e-7, 1e-9, 1e-6, big_m=False)


@milp
def test_asking_for_indicators_without_them_is_an_error(cycle_model):
    """Test that a solver without indicator constraints says so."""
    cycle_model.solver = "glpk"
    working, _ = prepare(cycle_model)

    assert not supports_indicators(working)
    with pytest.raises(ThermOptError, match="no native indicator"):
        solve_indicator_problem(
            working,
            tol=1e-4,
            tic_ids=[],
            basis=np.zeros((0, 0)),
            objective_ids={"R4"},
            maximise=True,
            use_indicators=True,
        )


@milp
def test_big_m_is_preferred_while_it_resolves_the_threshold(cycle_model):
    """Test that the cheaper formulation is chosen when it is adequate.

    The two agree on the answer and big-M is usually much the faster, so
    indicator constraints should only be reached for when the threshold is
    finer than big-M can resolve.

    """
    working, _ = prepare(cycle_model)
    working.tolerance = 1e-9

    # Bounds of 1000 and an integrality tolerance of 1e-9 leave 1e-6 of
    # slack, so the floor sits at SAFETY_FACTOR times that.
    assert _big_m_resolves(working, 1e-3, 1e-9, 1e-6)
    assert not _big_m_resolves(working, 1e-7, 1e-9, 1e-6)


@milp
def test_the_finer_threshold_switches_to_indicators(cycle_model, caplog):
    """Test that a threshold big-M cannot resolve is not simply rejected."""
    import logging

    working, _ = prepare(cycle_model)
    working.tolerance = 1e-9
    tic_ids, basis = cycle_basis(working, tol=1e-4, seed=0)

    with caplog.at_level(logging.INFO, logger="thermoptcobrapy._milp"):
        problem = solve_indicator_problem(
            working,
            tol=1e-6,
            tic_ids=tic_ids,
            basis=basis,
            objective_ids={rxn.id for rxn in working.reactions},
            maximise=True,
            minimum_for={rxn.id for rxn in working.reactions},
            time_limit=60,
        )

    assert "native indicator constraints" in caplog.text
    # and the answer is still right: the loop cannot carry flux
    assert set(problem.active_ids()) == {"R4", "EX_A", "EX_D"}


def test_a_fine_threshold_without_indicators_is_still_rejected(cycle_model):
    """Test that the fallback does not paper over an impossible request.

    A solver with no indicator constraints cannot resolve a threshold below
    the big-M floor, and should say so rather than return nonsense.

    """
    cycle_model.solver = "glpk"
    working, _ = prepare(cycle_model)

    with pytest.raises(ThermOptError, match="too small for this model"):
        solve_indicator_problem(
            working,
            tol=1e-7,
            tic_ids=[],
            basis=np.zeros((0, 0)),
            objective_ids={"R4"},
            maximise=True,
            minimum_for={"R4"},
        )

"""Fixtures for the ThermOpt tests.

The toy networks here are built so that the right answer is known by
inspection rather than by running the code being tested.
"""

import pytest
from cobra import Metabolite, Model, Reaction

#: The mixed-integer routines need a solver that actually solves MILPs.
MILP_SOLVER = "gurobi"


def _solver_available(name: str) -> bool:
    """Return whether cobrapy can use the named solver."""
    from cobra.util.solver import solvers

    return name in solvers


milp = pytest.mark.skipif(
    not _solver_available(MILP_SOLVER),
    reason=f"{MILP_SOLVER} is needed; GLPK does not solve these MILPs correctly",
)


def _build(reactions, exchanges, name="toy"):
    """Assemble a model from (id, stoichiometry) pairs."""
    model = Model(name)
    metabolites = {}

    def met(label):
        if label not in metabolites:
            metabolites[label] = Metabolite(label)
        return metabolites[label]

    built = []
    for rxn_id, stoichiometry in reactions:
        reaction = Reaction(rxn_id, lower_bound=-1000.0, upper_bound=1000.0)
        reaction.add_metabolites({met(k): v for k, v in stoichiometry.items()})
        built.append(reaction)
    for rxn_id, label in exchanges:
        reaction = Reaction(rxn_id, lower_bound=-1000.0, upper_bound=1000.0)
        reaction.add_metabolites({met(label): -1.0})
        built.append(reaction)
    model.add_reactions(built)
    model.solver = MILP_SOLVER
    return model


@pytest.fixture
def cycle_model() -> Model:
    """A->B->C->A, reachable only by cycling, plus a reaction that is not.

    ``R1``, ``R2`` and ``R3`` form a closed loop: ``B`` and ``C`` are touched
    by nothing else, so at steady state the three must carry identical flux
    and that flux is a cycle. Mass balance alone permits it; thermodynamics
    does not. ``R4`` drains ``A`` to ``D`` and is genuinely usable.

    """
    return _build(
        [
            ("R1", {"A": -1, "B": 1}),
            ("R2", {"B": -1, "C": 1}),
            ("R3", {"C": -1, "A": 1}),
            ("R4", {"A": -1, "D": 1}),
        ],
        [("EX_A", "A"), ("EX_D", "D")],
    )


@pytest.fixture
def double_cycle_model() -> Model:
    """The same network with a second loop sharing a reaction.

    Adding ``R5`` (B->A) makes three cycles: ``R1+R5``, ``R1+R2+R3``, and
    ``R2+R3`` against ``R5`` reversed.

    """
    return _build(
        [
            ("R1", {"A": -1, "B": 1}),
            ("R2", {"B": -1, "C": 1}),
            ("R3", {"C": -1, "A": 1}),
            ("R4", {"A": -1, "D": 1}),
            ("R5", {"B": -1, "A": 1}),
        ],
        [("EX_A", "A"), ("EX_D", "D")],
    )


@pytest.fixture
def acyclic_model() -> Model:
    """A linear path with no cycle at all."""
    return _build(
        [("R1", {"A": -1, "B": 1}), ("R2", {"B": -1, "C": 1})],
        [("EX_A", "A"), ("EX_C", "C")],
    )


@pytest.fixture
def reverse_only_model() -> Model:
    """A model whose only internal reaction is restricted to negative flux."""
    model = _build([("R1", {"A": -1, "B": 1})], [("EX_A", "A"), ("EX_B", "B")])
    model.reactions.R1.bounds = (-1000.0, 0.0)
    return model

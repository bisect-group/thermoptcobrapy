"""The thermodynamically constrained consistency check, ThermOptCC.

Flux consistency asks whether a reaction can carry flux at steady state.
That is a weaker question than it looks, because a reaction can be flux
consistent only with support from a TIC. ThermOptCC asks the stronger
question: can the reaction carry flux in a distribution that is also
thermodynamically feasible?

Reactions are categorised into four groups -- Forward, Reverse, Reversible
and Blocked -- because a reversible reaction may be thermodynamically
feasible in one direction and not the other.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, Optional, Tuple

from ._milp import solve_indicator_problem
from ._orientation import prepare
from ._tic import cycle_basis
from .exceptions import ThermOptError, ThermOptTimeLimitError

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)

BLOCKED = "Blocked"
UNBLOCKED = "Unblocked"
FORWARD = "Forward"
REVERSE = "Reverse"
REVERSIBLE = "Reversible"

#: Seconds allowed for the first solve.
INITIAL_TIME_LIMIT = 120.0
#: Added to the limit each time a solve runs out of time with nothing found.
TIME_LIMIT_STEP = 60.0
#: The limit stops growing here, so an unsolvable model raises rather than
#: running indefinitely.
MAX_TIME_LIMIT = 3600.0


def thermopt_cc(
    model: "Model",
    tol: float = 1e-4,
    get_direction: bool = False,
    seed: Optional[int] = None,
    time_limit: float = INITIAL_TIME_LIMIT,
    max_time_limit: float = MAX_TIME_LIMIT,
) -> Tuple[Dict[str, str], "Model"]:
    """Find the thermodynamically feasible flux directions of every reaction.

    Parameters
    ----------
    model : cobra.Model
        The model to check. It is left unchanged.
    tol : float, optional
        The flux a reaction has to reach to count as carrying flux
        (default 1e-4).
    get_direction : bool, optional
        Whether to report which directions are feasible rather than only
        whether the reaction is blocked (default False). Reporting directions
        costs roughly twice the work, because the reversible reactions are
        split and each half is tested separately.
    seed : int, optional
        A seed for the randomisation in both the consistency pre-check
        (SPECTRA's own) and Fast-SNP's random weights (default None).
    time_limit : float, optional
        Seconds allowed for the first mixed-integer solve of each round
        (default 120). Every round settles at least one reaction or this
        limit grows; it never shrinks back down within one call.
    max_time_limit : float, optional
        The ceiling `time_limit` grows to, by `TIME_LIMIT_STEP` seconds each
        time a solve times out having settled nothing new (default 3600).
        The MATLAB grows it without limit; this stops instead and raises, so
        an unsolvable model fails loudly rather than hanging or being
        reported as a result it was never proven to be.

    Returns
    -------
    tuple of (dict of {str: str}, cobra.Model)
        The verdict per reaction of the *returned* model, and that model: a
        copy in which the reactions irreversible in the reverse direction
        have been converted to forward irreversible reactions and the extreme
        bounds clamped. Verdicts are ``"Blocked"``/``"Unblocked"``, or with
        `get_direction` one of ``"Blocked"``, ``"Forward"``, ``"Reverse"``,
        ``"Reversible"``.

    Raises
    ------
    ThermOptError
        If the consistency program is infeasible outright -- typically a
        forced lower bound whose stoichiometry leaves no reaction room to be
        either exactly zero or at least `tol`.
    ThermOptTimeLimitError
        If a solve is still timing out once `time_limit` has grown to
        `max_time_limit`, with reactions left unsettled. Raise
        `max_time_limit`, or loosen `tol` so each solve is cheaper.

    Notes
    -----
    Results refer to the returned model, not to the one passed in, because
    the orientation step changes what "forward" means for any reaction that
    was irreversible in the reverse direction. This is the MATLAB's
    ``modModel``.

    """
    working, flipped = prepare(model)
    if flipped:
        logger.info("oriented %d reverse-only reactions forward", len(flipped))

    if get_direction:
        from ._irreversible import convert_to_irreversible

        searched, mapping = convert_to_irreversible(working)
    else:
        searched = working
        mapping = {rxn.id: [rxn.id] for rxn in working.reactions}

    tic_ids, basis = cycle_basis(searched, tol, seed=seed)
    logger.info(
        "%d reactions can cycle; Fast-SNP returned %d basis vectors",
        len(tic_ids),
        basis.shape[1] if basis.size else 0,
    )

    core = {rxn.id for rxn in searched.reactions}
    consistent = set()
    current_limit = time_limit
    first_pass = True

    while core:
        problem = solve_indicator_problem(
            searched,
            tol=tol,
            tic_ids=tic_ids,
            basis=basis,
            objective_ids=core,
            maximise=True,
            minimum_for=core,
            time_limit=current_limit,
        )
        if first_pass and problem.status == "infeasible":
            # Nothing is blocked here: no flux distribution at all satisfies
            # the program, so reporting every reaction as blocked would
            # present an infeasibility as a result.
            raise ThermOptError(
                f"The consistency program is infeasible at tol={tol:.3g}, so "
                f"no verdict can be given for any reaction. The usual cause "
                f"is a reaction with a forced lower bound -- a biomass "
                f"reaction pinned to a growth rate, say -- whose stoichiometry "
                f"demands flux below tol somewhere else: every reaction has to "
                f"be either exactly zero or at least tol, and that leaves no "
                f"room. Lower tol, or release the forced bounds if the "
                f"question is which reactions can carry flux at all."
            )
        first_pass = False

        found = {rxn_id for rxn_id in problem.active_ids() if rxn_id in core}
        if found:
            consistent |= found
            core -= found
            logger.debug(
                "%d newly consistent, %d still to settle", len(found), len(core)
            )
            continue

        if problem.status == "time_limit":
            if current_limit < max_time_limit:
                current_limit = min(current_limit + TIME_LIMIT_STEP, max_time_limit)
                logger.info(
                    "solve timed out; raising the limit to %gs", current_limit
                )
                continue
            raise ThermOptTimeLimitError(
                f"the consistency solve hit the {max_time_limit:g}s time limit "
                f"ceiling with {len(core)} reaction(s) still unsettled. Their "
                f"blocked/unblocked verdict was never proved, so reporting "
                f"them either way would be unsupported by the solve. Raise "
                f"max_time_limit, or loosen tol so each solve is cheaper."
            )
        break

    verdicts = {}
    for rxn in working.reactions:
        halves = mapping[rxn.id]
        live = [half in consistent for half in halves]
        if not get_direction:
            verdicts[rxn.id] = UNBLOCKED if any(live) else BLOCKED
        elif all(live) and len(live) == 2:
            verdicts[rxn.id] = REVERSIBLE
        elif not any(live):
            verdicts[rxn.id] = BLOCKED
        elif live[0]:
            verdicts[rxn.id] = FORWARD
        else:
            verdicts[rxn.id] = REVERSE

    return verdicts, working


def thermodynamically_blocked_ids(
    model: "Model", tol: float = 1e-4, seed: Optional[int] = None, **kwargs
) -> list:
    """Return the reactions that cannot carry thermodynamically feasible flux.

    Parameters
    ----------
    model : cobra.Model
        The model to check. It is left unchanged.
    tol : float, optional
        The flux threshold (default 1e-4).
    seed : int, optional
        A seed for Fast-SNP's random weights (default None).
    **kwargs
        Passed through to :func:`thermopt_cc`.

    Returns
    -------
    list of str
        The identifiers of the blocked reactions, in model order.

    """
    verdicts, working = thermopt_cc(model, tol=tol, seed=seed, **kwargs)
    return [rxn.id for rxn in working.reactions if verdicts[rxn.id] == BLOCKED]

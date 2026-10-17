"""Thermodynamically feasible context-specific models, ThermOptiCS.

Context-specific extraction normally asks for the smallest sub-network that
keeps a given set of core reactions active. That leaves a gap: "active" is
usually read as flux consistent, so a core reaction can be kept alive by a
cycle, and the model handed back inherits the cycle along with it.

ThermOptiCS closes the gap by asking for activity that is thermodynamically
feasible, using the same indicator program as the consistency check. It works
in rounds. Each round first establishes which of the still-unplaced core
reactions can be active together, then finds the cheapest set of non-core
reactions that supports them. Everything settled drops out of the core and
out of the penalty set, and the next round takes what is left.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Iterable, List, Optional, Tuple
from warnings import warn

from ._milp import solve_indicator_problem
from ._orientation import prepare
from ._subset import subset_model
from ._tic import cycle_basis
from .consistency import BLOCKED, thermopt_cc
from .exceptions import ThermOptInfeasibleCoreError

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)

#: Seconds allowed for each mixed-integer solve.
DEFAULT_TIME_LIMIT = 120.0


def thermoptics(
    model: "Model",
    core: Iterable[str],
    tol: float = 1e-4,
    do_tocc: bool = True,
    time_limit: float = DEFAULT_TIME_LIMIT,
    seed: Optional[int] = None,
) -> Tuple["Model", List[str]]:
    """Extract a thermodynamically feasible context-specific model.

    Parameters
    ----------
    model : cobra.Model
        The reference model to extract from. It is left unchanged.
    core : iterable of str
        The reactions the result has to keep active.
    tol : float, optional
        The flux a reaction has to reach to count as active (default 1e-4).
    do_tocc : bool, optional
        Whether to run the consistency check first and drop the
        thermodynamically blocked reactions from the core (default True).
        Pass False only when the model is already known to be
        thermodynamically flux consistent; a blocked core reaction then makes
        the extraction infeasible rather than merely slower.
    time_limit : float, optional
        Seconds allowed for each mixed-integer solve (default 120).
    seed : int, optional
        A seed for Fast-SNP's random weights (default None).

    Returns
    -------
    tuple of (cobra.Model, list of str)
        The extracted model, and the core reactions that had to be dropped
        because they are thermodynamically blocked.

    Raises
    ------
    ThermOptInfeasibleCoreError
        If no thermodynamically feasible flux distribution activates any of
        the remaining core reactions. With `do_tocc` left on this should not
        happen, since the blocked ones have already been removed. Unlike
        :func:`thermopt_cc`, a solve that merely times out before finding an
        incumbent is reported this way too, with the solver status in the
        message -- raising `time_limit` is the fix in that case, not
        `do_tocc`.

    Notes
    -----
    The result is built from the oriented working copy, so a reaction that
    was restricted to negative flux in the input appears with its
    stoichiometry negated and its bounds swapped. Its flux maps one to one
    onto the original's, negated.

    """
    working, flipped = prepare(model)
    if flipped:
        logger.info("oriented %d reverse-only reactions forward", len(flipped))

    all_ids = [rxn.id for rxn in working.reactions]
    known = set(all_ids)
    core = {rxn_id for rxn_id in core}
    unknown = core - known
    if unknown:
        raise KeyError(
            f"{len(unknown)} core reaction(s) are not in the model, "
            f"for example {sorted(unknown)[:5]}"
        )

    blocked_core: List[str] = []
    if do_tocc:
        verdicts, _ = thermopt_cc(
            working,
            tol=tol,
            get_direction=False,
            seed=seed,
            time_limit=time_limit,
        )
        blocked = {rxn_id for rxn_id, verdict in verdicts.items() if verdict == BLOCKED}
        blocked_core = [
            rxn_id for rxn_id in all_ids if rxn_id in core and rxn_id in blocked
        ]
        if blocked_core:
            warn(
                f"{len(blocked_core)} core reaction(s) are thermodynamically "
                "blocked and have been dropped from the core; they are "
                "returned alongside the model.",
                UserWarning,
                stacklevel=2,
            )
            core -= set(blocked_core)

    tic_ids, basis = cycle_basis(working, tol, seed=seed)
    logger.info(
        "%d reactions can cycle; Fast-SNP returned %d basis vectors",
        len(tic_ids),
        basis.shape[1] if basis.size else 0,
    )

    minimal: set = set()
    penalty = known - core
    remaining = set(core)

    while remaining:
        feasible = solve_indicator_problem(
            working,
            tol=tol,
            tic_ids=tic_ids,
            basis=basis,
            objective_ids=remaining,
            maximise=True,
            minimum_for=None,
            time_limit=time_limit,
        )
        active = {rxn_id for rxn_id in feasible.active_ids() if rxn_id in remaining}
        if not active:
            raise ThermOptInfeasibleCoreError(
                f"None of the {len(remaining)} remaining core reaction(s) can "
                f"carry thermodynamically feasible flux (solver status "
                f"{feasible.status!r}). Run with do_tocc=True, which drops "
                f"the blocked ones first, or raise time_limit."
            )

        smallest = solve_indicator_problem(
            working,
            tol=tol,
            tic_ids=tic_ids,
            basis=basis,
            objective_ids=penalty,
            maximise=False,
            minimum_for=None,
            force_active=active,
            time_limit=time_limit,
        )
        if not smallest.usable:
            raise ThermOptInfeasibleCoreError(
                f"The minimal-network step found nothing for {len(active)} "
                f"core reaction(s) already shown feasible (solver status "
                f"{smallest.status!r}). Raising time_limit usually helps."
            )

        # A reaction belongs to the network when its indicator is set or it
        # simply carries flux, which catches the support reactions the
        # objective had no reason to switch on explicitly.
        carried = {
            rxn_id
            for rxn_id, value in smallest.fluxes.items()
            if abs(value) >= tol * 1e-3
        }
        chosen = set(smallest.active_ids()) | carried
        minimal |= chosen
        remaining -= chosen
        penalty -= chosen
        logger.debug(
            "round kept %d reactions, %d core still to place",
            len(chosen),
            len(remaining),
        )

    minimal |= core & known
    extracted = subset_model(working, minimal)
    logger.info("extracted %d of %d reactions", len(extracted.reactions), len(all_ids))
    return extracted, blocked_core

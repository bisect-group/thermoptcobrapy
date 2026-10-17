"""Enumerating a model's thermodynamically infeasible cycles.

``thermopt_enum_milp`` minimises the *number* of active reactions in a cycle
sub-model and excludes each cycle it finds with a no-good cut, so the next
solve has to produce a different one. It needs a real mixed-integer solver.

It works on the cycle sub-model: exchanges shut, positive lower bounds
released, so the only flux left is cycling. It runs at a flux threshold of 1,
since in a closed network the scale of a cycle is arbitrary and only its
support matters.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Set, Tuple

from optlang.interface import FEASIBLE, OPTIMAL
from optlang.symbolics import Zero

from ._orientation import flip_reverse_only, normalise_bounds
from ._subset import subset_model
from ._tic import consistent_reaction_ids, exchange_reaction_ids

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)

#: The flux threshold used throughout the enumerators. A cycle can be scaled
#: freely, so this only sets the scale at which one is looked for.
CYCLE_TOL = 1.0
#: Fluxes below ``CYCLE_TOL * SUPPORT_FRACTION`` are not part of the cycle.
SUPPORT_FRACTION = 1e-3
USABLE = frozenset({OPTIMAL, FEASIBLE})


def _prepare_for_enumeration(model: "Model") -> "Model":
    """Return the closed, oriented copy both enumerators start from.

    Parameters
    ----------
    model : cobra.Model
        The model to prepare. It is left unchanged.

    Returns
    -------
    cobra.Model
        A copy with the reverse-only reactions oriented forward, every
        positive lower bound released and the extreme bounds clamped. The
        exchanges are left open: this is the model the returned cycles refer
        to, and it is what the MATLAB keeps as ``modModel``. Closing the
        exchanges is a separate step, done to the copy the search runs on.

    """
    working = model.copy()
    flip_reverse_only(working)
    for reaction in working.reactions:
        lower, upper = reaction.bounds
        if lower > 0.0:
            reaction.bounds = (0.0, upper)
    normalise_bounds(working)
    return working


def _close_exchanges(model: "Model") -> "Model":
    """Return a copy with every boundary reaction pinned shut.

    Parameters
    ----------
    model : cobra.Model
        The model to close. It is left unchanged.

    Returns
    -------
    cobra.Model
        A copy in which nothing enters or leaves, so the only flux it can
        carry is cycling.

    """
    closed = model.copy()
    exchanges = exchange_reaction_ids(closed)
    for reaction in closed.reactions:
        if reaction.id in exchanges:
            reaction.bounds = (0.0, 0.0)
    return closed


def reaction_order(model: "Model") -> List[str]:
    """Order reactions by the size of the component they belong to.

    Parameters
    ----------
    model : cobra.Model
        The model to order.

    Returns
    -------
    list of str
        The reaction identifiers, smallest connected component first and
        model order within a component.

    Notes
    -----
    Two reactions are connected when they share a metabolite. Taking the
    small components first means the enumerators finish whole components
    early and shrink the problem, rather than circling the largest one
    throughout.

    """
    neighbours: Dict[str, Set[str]] = {rxn.id: set() for rxn in model.reactions}
    for metabolite in model.metabolites:
        touching = [rxn.id for rxn in metabolite.reactions]
        for rxn_id in touching:
            neighbours[rxn_id].update(touching)

    seen: Set[str] = set()
    components: List[List[str]] = []
    for reaction in model.reactions:
        if reaction.id in seen:
            continue
        stack, component = [reaction.id], []
        seen.add(reaction.id)
        while stack:
            current = stack.pop()
            component.append(current)
            for other in neighbours[current]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        components.append(component)

    position = {rxn.id: i for i, rxn in enumerate(model.reactions)}
    components.sort(key=lambda comp: (len(comp), min(position[r] for r in comp)))
    return [
        rxn_id
        for component in components
        for rxn_id in sorted(component, key=lambda r: position[r])
    ]


def _force(reaction, direction: float, magnitude: float) -> bool:
    """Pin a reaction to carry flux in one direction, if it can.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to force.
    direction : float
        1 to force positive flux, -1 to force negative.
    magnitude : float
        The flux it has to reach.

    Returns
    -------
    bool
        False when the reaction's own bounds forbid that direction, in which
        case nothing was changed and the caller should treat the problem as
        infeasible rather than build it.

    """
    lower, upper = reaction.bounds
    if direction > 0:
        wanted = max(magnitude, lower)
        if wanted > upper:
            return False
        reaction.bounds = (wanted, upper)
    else:
        wanted = -magnitude
        if wanted < lower:
            return False
        reaction.bounds = (lower, wanted)
    return True


def _support(fluxes: Dict[str, float], order: Sequence[str]) -> List[str]:
    """Return the reactions carrying flux, in the given order.

    Parameters
    ----------
    fluxes : dict of {str: float}
        The flux of every reaction.
    order : sequence of str
        The order to report them in.

    Returns
    -------
    list of str
        The identifiers whose absolute flux reaches the support threshold.

    """
    cutoff = CYCLE_TOL * SUPPORT_FRACTION
    return [rxn_id for rxn_id in order if abs(fluxes.get(rxn_id, 0.0)) >= cutoff]


def _cycle_milp(
    model: "Model",
    core_id: str,
    direction: float,
    cuts: Sequence[Set[str]],
    minimise: bool,
    time_limit: Optional[float] = None,
) -> Optional[Dict[str, float]]:
    """Find a cycle through one reaction, excluding the cycles already found.

    Parameters
    ----------
    model : cobra.Model
        The cycle sub-model. Left unchanged.
    core_id : str
        The reaction the cycle has to pass through.
    direction : float
        1 or -1, the direction it is forced in.
    cuts : sequence of set of str
        Cycles already found. Each becomes a no-good cut forbidding all of
        its reactions from being active together, so the solver has to
        produce something new.
    minimise : bool
        Whether to minimise the number of active reactions. False asks only
        for feasibility, which is much cheaper and is used to find out
        whether any cycle is left at all.
    time_limit : float, optional
        Seconds allowed for the solve (default None).

    Returns
    -------
    dict of {str: float} or None
        The fluxes, or None when no cycle is left.

    """
    prob = model.problem
    indicators: Dict[str, object] = {}
    to_add, coefficients = [], []

    for reaction in model.reactions:
        if reaction.id == core_id:
            continue
        lower, upper = reaction.bounds
        switch = prob.Variable(f"to_y_{reaction.id}", type="binary")
        indicators[reaction.id] = switch
        to_add.append(switch)
        # v_i >= lb_i y_i and v_i <= ub_i y_i, so y_i off pins the reaction.
        low = prob.Constraint(Zero, name=f"to_yl_{reaction.id}", lb=0.0)
        high = prob.Constraint(Zero, name=f"to_yu_{reaction.id}", ub=0.0)
        to_add.extend([low, high])
        coefficients.append(
            (
                low,
                {
                    reaction.forward_variable: 1.0,
                    reaction.reverse_variable: -1.0,
                    switch: -lower,
                },
            )
        )
        coefficients.append(
            (
                high,
                {
                    reaction.forward_variable: 1.0,
                    reaction.reverse_variable: -1.0,
                    switch: -upper,
                },
            )
        )

    for number, cut in enumerate(cuts):
        terms = {indicators[r]: 1.0 for r in cut if r in indicators}
        if not terms:
            continue
        constraint = prob.Constraint(
            Zero, name=f"to_cut_{number}", ub=float(len(terms) - 1)
        )
        to_add.append(constraint)
        coefficients.append((constraint, terms))

    with model:
        if not _force(model.reactions.get_by_id(core_id), direction, CYCLE_TOL):
            return None
        model.add_cons_vars(to_add)
        model.solver.update()
        for constraint, terms in coefficients:
            constraint.set_linear_coefficients(terms)
        to_add.reverse()

        model.objective = prob.Objective(Zero, direction="min")
        if minimise and indicators:
            model.objective.set_linear_coefficients(
                {switch: 1.0 for switch in indicators.values()}
            )
        configuration = model.solver.configuration
        previous = configuration.timeout
        if time_limit is not None:
            configuration.timeout = int(max(1, round(time_limit)))
        try:
            model.solver.optimize()
            if model.solver.status not in USABLE:
                return None
            return {rxn.id: rxn.flux for rxn in model.reactions}
        finally:
            configuration.timeout = previous


def _tic_models(
    model: "Model",
    core_id: str,
    direction: float,
    cuts: List[Set[str]],
    time_limit: Optional[float] = None,
) -> Tuple[List[List[str]], List[List[float]]]:
    """Recursively pull out the cycles through one reaction.

    Parameters
    ----------
    model : cobra.Model
        The cycle sub-model. Left unchanged.
    core_id : str
        The reaction every cycle has to pass through.
    direction : float
        1 or -1.
    cuts : list of set of str
        Cycles already excluded. Extended as new ones are found.
    time_limit : float, optional
        Seconds per solve (default None).

    Returns
    -------
    tuple of (list of list of str, list of list of float)
        The cycles and their fluxes.

    Notes
    -----
    Each round first asks whether any cycle survives the cuts, then asks for
    the smallest one, then narrows the model to that cycle's reactions and
    recurses. Narrowing is what keeps the minimisation cheap: the next solve
    only sees a handful of reactions.

    """
    order = [rxn.id for rxn in model.reactions]
    feasible = _cycle_milp(
        model, core_id, direction, cuts, minimise=False, time_limit=time_limit
    )
    if feasible is None:
        return [], []
    support = _support(feasible, order)
    if core_id not in support:
        return [], []

    narrowed = subset_model(model, set(support))
    inner_cuts = [cut for cut in cuts if cut <= set(support)]
    smallest = _cycle_milp(
        narrowed,
        core_id,
        direction,
        inner_cuts,
        minimise=True,
        time_limit=time_limit,
    )
    if smallest is None:
        return [], []
    cycle = _support(smallest, [rxn.id for rxn in narrowed.reactions])
    if not cycle:
        return [], []

    fluxes = [smallest[rxn_id] for rxn_id in cycle]
    cuts.append(set(cycle))
    deeper, deeper_flux = _tic_models(
        narrowed, core_id, direction, cuts, time_limit=time_limit
    )
    return [cycle] + deeper, [fluxes] + deeper_flux


def thermopt_enum_milp(
    model: "Model",
    tol: float = CYCLE_TOL,
    time_limit: Optional[float] = None,
) -> Tuple[List[List[str]], List[List[float]], List[str], "Model"]:
    """Enumerate the infeasible cycles of a model by mixed-integer programming.

    Parameters
    ----------
    model : cobra.Model
        The model to search. It is left unchanged.
    tol : float, optional
        The flux threshold for the internal consistency check (default 1).
    time_limit : float, optional
        Seconds allowed for each solve (default None).

    Returns
    -------
    tuple of (list of list of str, list of list of float, list of str, cobra.Model)
        The cycles, the relative flux through each of their reactions
        normalised so the smallest is 1, every reaction taking part in a
        cycle, and the oriented working model the rest refer to.

    Notes
    -----
    Minimising the number of active reactions gives the cycles themselves,
    rather than a cycle plus whatever else the solver happened to leave on.
    Requires a real mixed-integer solver: Gurobi or CPLEX.

    """
    working = _prepare_for_enumeration(model)
    closed = _close_exchanges(working)
    cycling = consistent_reaction_ids(closed, tol)
    if not cycling:
        return [], [], [], working

    cycle_rxns = [rxn.id for rxn in working.reactions if rxn.id in cycling]
    search = subset_model(closed, cycling)
    pending = reaction_order(search)

    cycles: List[List[str]] = []
    directions: List[List[float]] = []

    while pending:
        head = pending[0]
        reversible = {rxn.id for rxn in search.reactions if rxn.lower_bound < 0.0}
        cuts: List[Set[str]] = []

        # Each pass narrows to the cycle it finds, so a pass can only reach
        # cycles inside that one. Repeating on the full model with the cuts
        # accumulated is what reaches the rest.
        forward_cycles: List[List[str]] = []
        while True:
            found, found_flux = _tic_models(
                search, head, 1.0, cuts, time_limit=time_limit
            )
            if not found:
                break
            for cycle, flux in zip(found, found_flux):
                cycles.append(cycle)
                directions.append(_normalise(flux))
                forward_cycles.append(cycle)

        # A cycle made entirely of reversible reactions runs backwards too,
        # so its mirror is recorded rather than searched for, and it is the
        # only thing the reverse pass has to be told to avoid.
        reverse_cuts: List[Set[str]] = []
        for cycle in forward_cycles:
            if set(cycle) <= reversible:
                position = cycles.index(cycle)
                cycles.append(list(cycle))
                directions.append([-value for value in directions[position]])
                reverse_cuts.append(set(cycle))

        while True:
            found, found_flux = _tic_models(
                search, head, -1.0, reverse_cuts, time_limit=time_limit
            )
            if not found:
                break
            for cycle, flux in zip(found, found_flux):
                cycles.append(cycle)
                directions.append(_normalise(flux))

        search.reactions.get_by_id(head).bounds = (0.0, 0.0)
        still = consistent_reaction_ids(search, tol)
        search = subset_model(search, still)
        pending = [rxn_id for rxn_id in pending if rxn_id in still]
        logger.debug("%d reactions left to enumerate", len(pending))

    logger.info("found %d cycles over %d reactions", len(cycles), len(cycle_rxns))
    return cycles, directions, cycle_rxns, working


def _normalise(fluxes: Sequence[float]) -> List[float]:
    """Scale a cycle's fluxes so the smallest magnitude is one.

    Parameters
    ----------
    fluxes : sequence of float
        The flux through each reaction of a cycle.

    Returns
    -------
    list of float
        The same proportions, scaled.

    """
    smallest = min(abs(value) for value in fluxes)
    if smallest == 0.0:
        return list(fluxes)
    return [value / smallest for value in fluxes]

"""The mixed-integer program the three ThermOpt searches share.

All three — the consistency check, the context-specific feasibility step and
the minimal-network step — solve the same problem with a different objective.
Its variables are

====================  ==========================================
:math:`v_i`           the flux, from the reaction's own bounds
:math:`z^+_i`         binary, "carries forward flux"
:math:`z^-_i`         binary, "carries reverse flux", reversible only
:math:`G_j`           continuous, a free-energy proxy per cycling reaction
====================  ==========================================

and its constraints, besides the model's own :math:`S v = 0`, are

.. math::

    v_i - ub_i z^+_i + \\underline{r}_i z^-_i &\\le 0 \\\\
    v_i - \\underline{f}_i z^+_i - lb_i z^-_i &\\ge 0 \\\\
    z^+_i + z^-_i &\\le 1 \\\\
    N_{snp}^\\mathsf{T} G &= 0 \\\\
    1001 z^+_k + G_j &\\le 1000 \\\\
    G_j - 1001 z^-_k &\\ge -1000

The first two are what make :math:`z` mean something: with both off the
reaction is pinned to zero flux, with :math:`z^+` on it must carry at least
:math:`\\underline{f}_i`, with :math:`z^-` on it must carry at most
:math:`-\\underline{r}_i`. The third stops a reaction claiming both directions
at once.

The last three are the thermodynamics, in the Fast-SNP form. :math:`G` plays
the role of a free-energy change: it is forced negative when a cycling
reaction runs forward and positive when it runs in reverse, and it is
confined to the orthogonal complement of the bounded null space. A cycle
would need :math:`G` to decrease all the way round and come back to itself,
which that complement does not allow, so no assignment of :math:`z` can
switch a whole cycle on.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from optlang.interface import FEASIBLE, OPTIMAL
from optlang.symbolics import Zero

from ._orientation import BOUND_NORM, is_reversible
from .exceptions import ThermOptError, ThermOptSolverError

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)

#: The big-M coefficient on the indicator in the free-energy constraints,
#: one more than the bound on G so that the constraint is slack when the
#: indicator is off.
BIG_M = BOUND_NORM + 1.0
#: The tightest feasibility tolerance the *big-M* programs may be solved at.
#:
#: The big-M coefficients reach the largest flux bound, so an indicator is
#: only worth anything while ``M * integrality`` stays well under ``tol``.
#: That pushes the integrality tolerance down, and the obvious next step --
#: tightening feasibility to match, as the MATLAB does with feasTol 1e-9 --
#: turns out to break Gurobi 13 on these problems: it returns "optimal" with
#: every indicator off while a solution with hundreds of them on is feasible
#: to 1e-12, verified outside the solver. At 1e-6 the same problems solve
#: correctly and monotonically.
#:
#: This is a property of the big-M form, not of the solver: with native
#: indicator constraints the same model at 1e-9 returns the same 1367
#: reactions it returns at 1e-6. So the floor applies only when big-M is in
#: use, and the native path is free to work to whatever tolerance the model
#: asks for.
MIN_FEASIBILITY = 1e-6

#: How far ``tol`` must clear the slack the solver is allowed on an
#: indicator before the indicator is worth anything.
#:
#: Measured, not assumed. On an E. coli model with bounds of 1000 and an
#: integrality tolerance of 1e-9 the slack is 1e-6; thresholds of 2e-5 and
#: 3e-5 (20x and 30x that) still produced solutions contradicting their own
#: indicators, while 5e-5 (50x) solved correctly. The factor sits between
#: them: a factor of ten, the obvious choice, is not enough margin here.
SAFETY_FACTOR = 40.0

#: How far a reported flux may drift from what its indicator says before the
#: audit calls it a contradiction. A flux is a difference of two variables,
#: each held only to the feasibility tolerance, so some drift is expected.
#: Kept well under ``SAFETY_FACTOR`` so the audit still bites.
AUDIT_FACTOR = 10.0

#: Statuses that come with a usable incumbent.
USABLE = frozenset({OPTIMAL, FEASIBLE, "time_limit", "aborted"})


class _Problem:
    """The assembled program, held so the caller can read its solution back.

    Attributes
    ----------
    status : str
        The solver status reached.
    fluxes : dict of {str: float}
        The flux of every reaction, empty when nothing was found.
    forward, reverse : dict of {str: float}
        The ``z+`` and ``z-`` values, empty when nothing was found.

    """

    def __init__(self, status, fluxes, forward, reverse):
        self.status = status
        self.fluxes = fluxes
        self.forward = forward
        self.reverse = reverse

    @property
    def usable(self) -> bool:
        """Whether a solution was found at all."""
        return bool(self.forward)

    def active_ids(self, threshold: float = 0.5) -> List[str]:
        """Return the reactions switched on in either direction.

        Parameters
        ----------
        threshold : float, optional
            The value above which a binary counts as set (default 0.5).

        Returns
        -------
        list of str
            The identifiers, in the order the model had them.

        """
        return [
            rxn_id
            for rxn_id in self.forward
            if self.forward[rxn_id] >= threshold
            or self.reverse.get(rxn_id, 0.0) >= threshold
        ]


def _big_m_resolves(
    model: "Model", tol: float, integrality: float, feasibility: float
) -> bool:
    """Return whether the big-M formulation can resolve this threshold.

    Parameters
    ----------
    model : cobra.Model
        The model whose bounds set the big-M coefficients.
    tol : float
        The flux threshold wanted.
    integrality : float
        The integrality tolerance the solve would use.
    feasibility : float
        The feasibility tolerance the solve would use.

    Returns
    -------
    bool
        False when :func:`check_tolerances` would reject the combination.

    """
    try:
        check_tolerances(model, tol, integrality, feasibility, big_m=True)
    except ThermOptError:
        return False
    return True


def supports_indicators(model: "Model") -> bool:
    """Return whether the model's solver offers native indicator constraints.

    Parameters
    ----------
    model : cobra.Model
        The model whose solver to inspect.

    Returns
    -------
    bool
        True for Gurobi. CPLEX has them too but is not wired up here, so it
        takes the big-M path until someone with a licence can test it.

    """
    return "gurobi" in type(model.solver).__module__.lower()


def _attach_indicators(
    model: "Model",
    reactions,
    forward_z,
    reverse_z,
    either_z,
    energy,
    tic_ids,
    minima,
):
    """Attach the indicator constraints to the solver, bypassing optlang.

    Parameters
    ----------
    model : cobra.Model
        The model to attach to. Its problem must already hold the variables.
    reactions : list of cobra.Reaction
        The reactions to constrain.
    forward_z, reverse_z, either_z : dict of {str: optlang.Variable}
        The indicator variables.
    energy : dict of {str: optlang.Variable}
        The free-energy variables, keyed by cycling reaction.
    tic_ids : sequence of str
        The cycling reactions.
    minima : dict of {str: tuple of (float, float)}
        The forward and reverse minima per reaction.

    Returns
    -------
    list
        Handles for the constraints added, to be removed afterwards.

    Notes
    -----
    An indicator constraint says "when this binary is set, that inequality
    holds", and the solver enforces it by branching rather than by arithmetic.
    There is no coefficient to multiply a not-quite-integral binary by, so the
    ``M * integrality`` slack that limits the big-M form simply does not
    arise, and ``tol`` is free to go as low as the feasibility tolerance
    allows.

    optlang has no portable way to express these, so they are added straight
    to the underlying solver and removed again before the surrounding context
    tears the variables down.

    """
    from gurobipy import GRB

    problem = model.solver.problem
    problem.update()

    def native(variable):
        found = problem.getVarByName(variable.name)
        if found is None:
            raise ThermOptError(
                f"the solver has no variable named {variable.name!r}; the "
                "problem was not updated before the indicators were attached"
            )
        return found

    added = []
    for reaction in reactions:
        rxn_id = reaction.id
        flux = native(reaction.forward_variable) - native(reaction.reverse_variable)
        min_fwd, min_rev = minima[rxn_id]
        zp = native(forward_z[rxn_id])
        added.append(
            problem.addGenConstrIndicator(
                zp, True, flux, GRB.GREATER_EQUAL, min_fwd, f"ti_f_{rxn_id}"
            )
        )
        if rxn_id in reverse_z:
            zm = native(reverse_z[rxn_id])
            added.append(
                problem.addGenConstrIndicator(
                    zm, True, flux, GRB.LESS_EQUAL, -min_rev, f"ti_r_{rxn_id}"
                )
            )
            # Neither direction set means no flux at all.
            added.append(
                problem.addGenConstrIndicator(
                    native(either_z[rxn_id]), False, flux, GRB.EQUAL, 0.0,
                    f"ti_z_{rxn_id}",
                )
            )
        else:
            added.append(
                problem.addGenConstrIndicator(
                    zp, False, flux, GRB.EQUAL, 0.0, f"ti_z_{rxn_id}"
                )
            )

    # The free energy falls when a cycling reaction runs forward and rises
    # when it runs backward, which is what forbids a closed loop.
    for rxn_id in tic_ids:
        potential = native(energy[rxn_id])
        added.append(
            problem.addGenConstrIndicator(
                native(forward_z[rxn_id]), True, potential, GRB.LESS_EQUAL,
                -1.0, f"ti_gf_{rxn_id}",
            )
        )
        if rxn_id in reverse_z:
            added.append(
                problem.addGenConstrIndicator(
                    native(reverse_z[rxn_id]), True, potential,
                    GRB.GREATER_EQUAL, 1.0, f"ti_gr_{rxn_id}",
                )
            )
    problem.update()
    return added


def _detach_indicators(model: "Model", handles) -> None:
    """Remove indicator constraints added outside optlang's bookkeeping.

    Parameters
    ----------
    model : cobra.Model
        The model they were attached to.
    handles : list
        What :func:`_attach_indicators` returned.

    Notes
    -----
    This has to happen before the surrounding context removes the variables
    the constraints refer to, or the solver is left holding references to
    things that no longer exist.

    """
    if not handles:
        return
    problem = model.solver.problem
    problem.remove(handles)
    problem.update()


def check_tolerances(
    model: "Model",
    tol: float,
    integrality: float,
    feasibility: float,
    big_m: bool = True,
) -> None:
    """Raise if the flux threshold is too small for the big-M formulation.

    Parameters
    ----------
    model : cobra.Model
        The model whose bounds set the big-M coefficients.
    tol : float
        The flux a reaction has to reach to count as carrying flux.
    integrality : float
        The integrality tolerance the solve will use.
    feasibility : float
        The feasibility tolerance the solve will use.

    Raises
    ------
    ThermOptError
        If ``tol`` is not comfortably above the slack the solver may take.

    Notes
    -----
    An indicator is enforced through a term like ``ub_i * z_i``, so a binary
    sitting ``integrality`` away from zero buys ``ub_i * integrality`` of
    flux for free. Once that reaches ``tol`` the indicator stops meaning
    anything and the solver returns answers that are wrong rather than
    imprecise -- a reaction marked active while carrying nothing. The same
    goes for a flux smaller than the feasibility tolerance, which the solver
    cannot distinguish from zero.

    Gurobi will not accept an integrality tolerance below 1e-9, so on a model
    whose bounds reach 1000 the smallest workable threshold is around 5e-5.
    Wanting to resolve finer flux than that means the big-M formulation is
    the wrong tool, not that the threshold should be forced down.

    The bound that matters is the largest one in the model, so narrowing the
    flux bounds of the handful of reactions left at a default 1000 lowers the
    floor for everything.

    """
    largest = max(
        (max(abs(rxn.lower_bound), abs(rxn.upper_bound)) for rxn in model.reactions),
        default=0.0,
    )
    phantom = largest * integrality if big_m else 0.0
    if tol < SAFETY_FACTOR * phantom:
        raise ThermOptError(
            f"tol={tol:.3g} is too small for this model: its largest flux "
            f"bound is {largest:g}, and with an integrality tolerance of "
            f"{integrality:.3g} the solver may grant {phantom:.3g} of flux to "
            f"a reaction whose indicator is off. Raise tol to at least "
            f"{SAFETY_FACTOR * phantom:.3g}, or narrow the model's flux "
            f"bounds. Reactions that cannot carry {tol:.3g} are reported "
            f"blocked, which is the correct result at this threshold."
        )
    if tol < SAFETY_FACTOR * feasibility:
        raise ThermOptError(
            f"tol={tol:.3g} is below {SAFETY_FACTOR:g}x the solver's "
            f"feasibility tolerance ({feasibility:.3g}), so a flux of tol "
            f"cannot be told apart from zero. Raise tol to at least "
            f"{SAFETY_FACTOR * feasibility:.3g}."
        )


def _audit(
    fluxes: Dict[str, float],
    forward: Dict[str, float],
    reverse: Dict[str, float],
    minima: Dict[str, Tuple[float, float]],
    atol: float,
) -> Optional[str]:
    """Return the first indicator the solution contradicts, or None.

    Parameters
    ----------
    fluxes : dict of {str: float}
        The reported flux of every reaction.
    forward, reverse : dict of {str: float}
        The reported indicator values.
    minima : dict of {str: tuple of (float, float)}
        The forward and reverse minima each reaction was held to.
    atol : float
        How far a value may stray before it counts as a contradiction.

    Returns
    -------
    str or None
        A description of the first violation found.

    Notes
    -----
    The indicators are the answer this package reports, so a solver that
    returns them inconsistent with the fluxes returns a wrong answer rather
    than an imprecise one. GLPK does exactly that on these problems: it has
    been seen reporting ``optimal`` with an indicator set on a reaction
    carrying no flux. Checking is cheap; trusting is not.

    """
    for rxn_id, flux in fluxes.items():
        min_fwd, min_rev = minima[rxn_id]
        on_forward = forward.get(rxn_id, 0.0) >= 0.5
        on_reverse = reverse.get(rxn_id, 0.0) >= 0.5
        if on_forward and flux < min_fwd - atol:
            return (
                f"{rxn_id} is marked forward but carries {flux:.3g}, "
                f"under the {min_fwd:.3g} required"
            )
        if on_reverse and flux > -min_rev + atol:
            return (
                f"{rxn_id} is marked reverse but carries {flux:.3g}, "
                f"over the {-min_rev:.3g} required"
            )
        if not on_forward and not on_reverse and abs(flux) > atol:
            return f"{rxn_id} is marked inactive but carries {flux:.3g}"
    return None


def solve_indicator_problem(
    model: "Model",
    tol: float,
    tic_ids: Sequence[str],
    basis: np.ndarray,
    objective_ids: Iterable[str],
    maximise: bool,
    minimum_for: Optional[Iterable[str]] = None,
    force_active: Optional[Iterable[str]] = None,
    time_limit: Optional[float] = None,
    integrality: float = 1e-9,
    use_indicators: Optional[bool] = None,
) -> _Problem:
    """Build and solve the indicator program on a model.

    Parameters
    ----------
    model : cobra.Model
        The model to operate on, already oriented. Left unchanged: everything
        is added inside a context and rolled back.
    tol : float
        The flux a reaction must reach to count as carrying flux.
    tic_ids : sequence of str
        The cycling reactions, in model order, matching the rows of `basis`.
    basis : numpy.ndarray
        The Fast-SNP basis, ``len(tic_ids) x K``. May be empty, in which case
        no free-energy variables are added and the program is a plain
        indicator model.
    objective_ids : iterable of str
        The reactions whose indicators carry an objective coefficient of one.
    maximise : bool
        Whether to maximise that sum rather than minimise it.
    minimum_for : iterable of str, optional
        The reactions that must reach `tol` when their indicator is set. The
        default, None, applies it to every reaction; passing a subset is what
        the consistency check does, where only the core needs forcing.
    force_active : iterable of str, optional
        Reactions required to carry flux in one direction or the other
        (default None).
    time_limit : float, optional
        Seconds allowed for the solve (default None, the solver's own).
    integrality : float, optional
        The integrality tolerance (default 1e-9, as the MATLAB sets).
    use_indicators : bool, optional
        Whether to express the indicators natively rather than through big-M
        coefficients. The default, None, picks whichever fits: big-M while
        `tol` clears the slack it leaves, native indicators when it does not
        and the solver offers them. Big-M is the faster of the two on most
        models, often by a wide margin, and the two return the same verdicts,
        so there is no reason to pay for indicators until the threshold
        demands them. Forcing either is allowed; forcing True on a solver
        without them is an error.

    Returns
    -------
    _Problem
        The solution, or an empty one when the solver found nothing.

    """
    own_tolerance = getattr(model, "tolerance", MIN_FEASIBILITY)

    if use_indicators is None:
        # Prefer big-M, which is usually much the faster, and reach for
        # indicator constraints only when the threshold is finer than big-M
        # can resolve. The two agree on the answer, so this costs nothing.
        native = supports_indicators(model) and not _big_m_resolves(
            model, tol, integrality, max(MIN_FEASIBILITY, own_tolerance)
        )
        if native:
            logger.info(
                "tol=%.3g is finer than big-M can resolve on this model; "
                "using native indicator constraints instead",
                tol,
            )
    else:
        native = use_indicators
        if native and not supports_indicators(model):
            raise ThermOptError(
                f"{type(model.solver).__module__} has no native indicator "
                "constraints; leave use_indicators unset to fall back to "
                "big-M."
            )

    # Only big-M needs protecting from a tight feasibility tolerance.
    feasibility = own_tolerance if native else max(MIN_FEASIBILITY, own_tolerance)
    check_tolerances(model, tol, integrality, feasibility, big_m=not native)

    prob = model.problem
    reactions = list(model.reactions)
    objective_ids = set(objective_ids)
    forced = set(force_active or ())
    needs_minimum = (
        set(minimum_for) if minimum_for is not None else {rxn.id for rxn in reactions}
    )

    n_basis = basis.shape[1] if basis.size else 0

    forward_z: Dict[str, object] = {}
    reverse_z: Dict[str, object] = {}
    either_z: Dict[str, object] = {}
    energy: Dict[str, object] = {}
    to_add = []
    coefficients: List[Tuple[object, Dict[object, float]]] = []

    for reaction in reactions:
        zp = prob.Variable(f"to_zp_{reaction.id}", type="binary")
        forward_z[reaction.id] = zp
        to_add.append(zp)
        if is_reversible(reaction):
            zm = prob.Variable(f"to_zm_{reaction.id}", type="binary")
            reverse_z[reaction.id] = zm
            to_add.append(zm)

    for rxn_id in tic_ids:
        variable = prob.Variable(f"to_g_{rxn_id}", lb=-BOUND_NORM, ub=BOUND_NORM)
        energy[rxn_id] = variable
        to_add.append(variable)

    minima: Dict[str, Tuple[float, float]] = {}
    for reaction in reactions:
        lower, upper = reaction.bounds
        minima[reaction.id] = (
            max(tol, lower) if reaction.id in needs_minimum else 0.0,
            tol if reaction.id in needs_minimum else 0.0,
        )

    # The nullspace equality is the same either way: it carries no indicator
    # and so needs no big-M.
    for column in range(n_basis):
        terms = {
            energy[tic_ids[row]]: float(basis[row, column])
            for row in range(len(tic_ids))
            if basis[row, column] != 0.0
        }
        if not terms:
            continue
        constraint = prob.Constraint(Zero, name=f"to_null_{column}", lb=0.0, ub=0.0)
        to_add.append(constraint)
        coefficients.append((constraint, terms))

    if native:
        # Only the pieces that cannot be expressed as an indicator are added
        # here; the rest is attached to the solver directly once the problem
        # exists, by _attach_indicators.
        for reaction in reactions:
            zm = reverse_z.get(reaction.id)
            if zm is None:
                continue
            # z_any = z+ + z-, so that "neither direction" is a single binary
            # an indicator can key off. It also subsumes z+ + z- <= 1.
            switch = prob.Variable(f"to_zany_{reaction.id}", type="binary")
            either_z[reaction.id] = switch
            to_add.append(switch)
            constraint = prob.Constraint(
                Zero, name=f"to_any_{reaction.id}", lb=0.0, ub=0.0
            )
            to_add.append(constraint)
            coefficients.append(
                (constraint, {switch: 1.0, forward_z[reaction.id]: -1.0, zm: -1.0})
            )
        for rxn_id in forced:
            terms = {forward_z[rxn_id]: 1.0}
            if rxn_id in reverse_z:
                terms[reverse_z[rxn_id]] = 1.0
            constraint = prob.Constraint(Zero, name=f"to_core_{rxn_id}", lb=1.0)
            to_add.append(constraint)
            coefficients.append((constraint, terms))
    else:
        for reaction in reactions:
            rxn_id = reaction.id
            zp = forward_z[rxn_id]
            zm = reverse_z.get(rxn_id)
            lower, upper = reaction.bounds
            min_fwd, min_rev = minima[rxn_id]

            # v_i - ub_i z+_i + min_rev_i z-_i <= 0
            upper_terms = {
                reaction.forward_variable: 1.0,
                reaction.reverse_variable: -1.0,
                zp: -upper,
            }
            if zm is not None:
                upper_terms[zm] = min_rev
            constraint = prob.Constraint(Zero, name=f"to_up_{rxn_id}", ub=0.0)
            to_add.append(constraint)
            coefficients.append((constraint, upper_terms))

            # v_i - min_fwd_i z+_i - lb_i z-_i >= 0
            lower_terms = {
                reaction.forward_variable: 1.0,
                reaction.reverse_variable: -1.0,
                zp: -min_fwd,
            }
            if zm is not None:
                lower_terms[zm] = -lower
            constraint = prob.Constraint(Zero, name=f"to_lo_{rxn_id}", lb=0.0)
            to_add.append(constraint)
            coefficients.append((constraint, lower_terms))

            # z+_i + z-_i <= 1
            if zm is not None:
                constraint = prob.Constraint(Zero, name=f"to_one_{rxn_id}", ub=1.0)
                to_add.append(constraint)
                coefficients.append((constraint, {zp: 1.0, zm: 1.0}))

            if rxn_id in forced:
                terms = {zp: 1.0}
                if zm is not None:
                    terms[zm] = 1.0
                constraint = prob.Constraint(Zero, name=f"to_core_{rxn_id}", lb=1.0)
                to_add.append(constraint)
                coefficients.append((constraint, terms))

        # The free energy must fall when a cycling reaction runs forward and
        # rise when it runs backward, which is what forbids a closed loop.
        for rxn_id in tic_ids:
            variable = energy[rxn_id]
            constraint = prob.Constraint(Zero, name=f"to_gf_{rxn_id}", ub=BOUND_NORM)
            to_add.append(constraint)
            coefficients.append(
                (constraint, {variable: 1.0, forward_z[rxn_id]: BIG_M})
            )
            zm = reverse_z.get(rxn_id)
            if zm is not None:
                constraint = prob.Constraint(
                    Zero, name=f"to_gr_{rxn_id}", lb=-BOUND_NORM
                )
                to_add.append(constraint)
                coefficients.append((constraint, {variable: 1.0, zm: -BIG_M}))

    status = None
    fluxes: Dict[str, float] = {}
    forward_values: Dict[str, float] = {}
    reverse_values: Dict[str, float] = {}

    with model:
        model.add_cons_vars(to_add)
        model.solver.update()
        for constraint, terms in coefficients:
            constraint.set_linear_coefficients(terms)
        # Reversed so the rollback deletes from the tail of optlang's
        # containers rather than re-indexing from the front.
        to_add.reverse()

        handles = []
        if native:
            handles = _attach_indicators(
                model,
                reactions,
                forward_z,
                reverse_z,
                either_z,
                energy,
                tic_ids,
                minima,
            )

        objective_terms = {}
        for rxn_id in objective_ids:
            if rxn_id in forward_z:
                objective_terms[forward_z[rxn_id]] = 1.0
            if rxn_id in reverse_z:
                objective_terms[reverse_z[rxn_id]] = 1.0
        model.objective = prob.Objective(Zero, direction="max" if maximise else "min")
        if objective_terms:
            model.objective.set_linear_coefficients(objective_terms)

        configuration = model.solver.configuration
        previous_timeout = configuration.timeout
        previous_integrality = previous_feasibility = None
        try:
            previous_integrality = configuration.tolerances.integrality
            configuration.tolerances.integrality = integrality
        except (AttributeError, ValueError):
            logger.debug("this solver does not expose an integrality tolerance")
        try:
            previous_feasibility = configuration.tolerances.feasibility
            if not native and previous_feasibility < MIN_FEASIBILITY:
                logger.info(
                    "raising the feasibility tolerance from %g to %g for this "
                    "mixed-integer solve; tighter than that, the solver "
                    "misreports the optimum",
                    previous_feasibility,
                    MIN_FEASIBILITY,
                )
                configuration.tolerances.feasibility = MIN_FEASIBILITY
        except (AttributeError, ValueError):
            logger.debug("this solver does not expose a feasibility tolerance")
        if time_limit is not None:
            # glpk's binding wants a whole number of seconds.
            configuration.timeout = int(max(1, round(time_limit)))
        try:
            model.solver.optimize()
            status = model.solver.status
            if status in USABLE:
                try:
                    forward_values = {
                        rxn_id: var.primal for rxn_id, var in forward_z.items()
                    }
                    reverse_values = {
                        rxn_id: var.primal for rxn_id, var in reverse_z.items()
                    }
                    fluxes = {rxn.id: rxn.flux for rxn in reactions}
                except Exception:  # noqa: BLE001 - no incumbent to read
                    forward_values, reverse_values, fluxes = {}, {}, {}
                if fluxes:
                    # A flux is the difference of two variables, each only
                    # held to the feasibility tolerance, so a few multiples
                    # of it is drift rather than a violation. check_tolerances
                    # has already guaranteed tol is at least this far above
                    # the noise, so the audit still bites on a real one.
                    atol = AUDIT_FACTOR * feasibility
                    complaint = _audit(
                        fluxes, forward_values, reverse_values, minima, atol
                    )
                    if complaint is not None:
                        raise ThermOptSolverError(
                            f"{type(model.solver).__module__} returned a "
                            f"solution its own indicator constraints forbid "
                            f"({complaint}). This is a solver defect, not a "
                            f"modelling one; GLPK is known to do it on these "
                            f"mixed-integer problems. Use Gurobi or CPLEX.",
                            status=status,
                        )
        finally:
            # Before the context removes the variables these refer to.
            _detach_indicators(model, handles)
            configuration.timeout = previous_timeout
            for attribute, value in (
                ("integrality", previous_integrality),
                ("feasibility", previous_feasibility),
            ):
                if value is None:
                    continue
                try:
                    setattr(configuration.tolerances, attribute, value)
                except (AttributeError, ValueError):
                    pass

    return _Problem(status, fluxes, forward_values, reverse_values)

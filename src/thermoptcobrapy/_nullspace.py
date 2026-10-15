r"""Fast Sparse Nullspace Pursuit.

A port of ``fast_snp.m``, itself Pedro A. Saa's Fast-SNP. It returns a sparse
basis for the part of the null space that the flux bounds actually allow,
which is what the loopless constraints are written against: a thermodynamically
infeasible cycle lives in the null space of ``S`` and respects the bounds, so a
basis of that subspace is enough to exclude every cycle at once.

Each basis vector comes from a pair of LPs,

.. math::

    \\min \; \\|v\\|_1 \\quad\\text{s.t.}\\quad S v = 0,\;
    x^\\mathsf{T} v \\ge \\varepsilon,\; lb \\le v \\le ub

and its mirror with :math:`-x`, where :math:`x` is a randomly weighted
projection away from everything found so far. Minimising the one-norm is what
makes the basis sparse; projecting is what keeps each vector independent of
its predecessors. The shorter of the two solutions is kept.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Optional

import numpy as np
from optlang.interface import OPTIMAL
from optlang.symbolics import Zero

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)

#: The right-hand side of the projection constraint. Large enough that a
#: solution cannot be a numerically scaled nothing, small enough not to clash
#: with the bounds.
EPSILON = 1e-1
#: The upper bound on the one-norm auxiliaries, as in ``buildSNVproblem``.
NORM_BOUND = 1e2


def _orthonormal_basis(matrix: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Return an orthonormal basis for the column space of a matrix.

    Parameters
    ----------
    matrix : numpy.ndarray
        The matrix whose column space to span.
    tol : float, optional
        Singular values at or below this count as zero (default 1e-12).

    Returns
    -------
    numpy.ndarray
        A matrix whose columns are an orthonormal basis, equivalent to
        MATLAB's ``orth``.

    """
    if matrix.size == 0:
        return np.zeros((matrix.shape[0], 0))
    left, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    cutoff = (
        max(matrix.shape)
        * np.finfo(float).eps
        * (singular[0] if singular.size else 0.0)
    )
    keep = singular > max(cutoff, tol)
    return left[:, keep]


def nullity(stoichiometry: np.ndarray) -> int:
    """Return the dimension of the null space of a stoichiometric matrix.

    Parameters
    ----------
    stoichiometry : numpy.ndarray
        The matrix to inspect.

    Returns
    -------
    int
        ``n - rank(S)``, the number of basis vectors Fast-SNP will look for.

    Notes
    -----
    Computed densely, as MATLAB's ``null(full(S))`` is. That is affordable
    here because Fast-SNP is only ever run on the cycle sub-model, which is a
    small fraction of a genome-scale network, never on the whole model.

    """
    if stoichiometry.size == 0:
        return 0
    return stoichiometry.shape[1] - int(np.linalg.matrix_rank(stoichiometry))


def fast_snp(
    model: "Model",
    v_tol: float = 1e-6,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Return a sparse basis of the bounded null space of a model.

    Parameters
    ----------
    model : cobra.Model
        The model to work on, normally the cycle sub-model. It is not
        modified.
    v_tol : float, optional
        Fluxes below this magnitude are read as zero when the sparsity of a
        candidate vector is judged (default 1e-6).
    seed : int, optional
        A seed for the random weights, making the basis reproducible
        (default None).

    Returns
    -------
    numpy.ndarray
        An ``n x K`` array whose rows follow ``model.reactions`` and whose
        columns are the basis vectors. ``K`` is at most the nullity of ``S``
        and is often smaller, since the search stops as soon as the bounds
        admit nothing new.

    """
    from cobra.util.array import create_stoichiometric_matrix

    stoichiometry = np.asarray(
        create_stoichiometric_matrix(model, array_type="dense"), dtype=float
    )
    n_reactions = stoichiometry.shape[1]
    if n_reactions == 0:
        return np.zeros((0, 0))

    target = nullity(stoichiometry)
    if target == 0:
        return np.zeros((n_reactions, 0))

    rng = np.random.default_rng(seed)
    prob = model.problem

    # One LP, reused. Everything but the projection row is fixed, so only its
    # coefficients are rewritten between iterations; rebuilding the whole
    # problem once per basis vector is the obvious way and much the slower.
    lp = prob.Model(name="fast_snp")
    flux_vars, norm_vars = [], []
    for reaction in model.reactions:
        flux_vars.append(
            prob.Variable(
                f"v_{reaction.id}",
                lb=reaction.lower_bound,
                ub=reaction.upper_bound,
            )
        )
        norm_vars.append(prob.Variable(f"t_{reaction.id}", lb=0.0, ub=NORM_BOUND))
    lp.add(flux_vars + norm_vars)
    lp.update()

    constraints = []
    for row, metabolite in enumerate(model.metabolites):
        terms = {
            flux_vars[col]: stoichiometry[row, col]
            for col in range(n_reactions)
            if stoichiometry[row, col] != 0.0
        }
        if terms:
            constraints.append(
                (
                    prob.Constraint(Zero, name=f"mb_{metabolite.id}", lb=0.0, ub=0.0),
                    terms,
                )
            )
    for index, reaction in enumerate(model.reactions):
        # t_i >= v_i and t_i >= -v_i, together t_i >= |v_i|.
        constraints.append(
            (
                prob.Constraint(Zero, name=f"nupper_{reaction.id}", lb=0.0),
                {norm_vars[index]: 1.0, flux_vars[index]: -1.0},
            )
        )
        constraints.append(
            (
                prob.Constraint(Zero, name=f"nlower_{reaction.id}", lb=0.0),
                {norm_vars[index]: 1.0, flux_vars[index]: 1.0},
            )
        )
    projection = prob.Constraint(Zero, name="projection", lb=EPSILON)
    constraints.append((projection, {}))

    lp.add([constraint for constraint, _ in constraints])
    lp.update()
    for constraint, terms in constraints:
        if terms:
            constraint.set_linear_coefficients(terms)

    lp.objective = prob.Objective(Zero, direction="min")
    lp.objective.set_linear_coefficients({var: 1.0 for var in norm_vars})

    basis = np.zeros((n_reactions, 0))
    for _ in range(target):
        weights = 1.0 + rng.random(n_reactions)
        if basis.shape[1]:
            orthonormal = _orthonormal_basis(basis)
            weights = weights - weights @ orthonormal @ orthonormal.T

        candidates = []
        for sign in (1.0, -1.0):
            projection.set_linear_coefficients(
                {flux_vars[i]: sign * weights[i] for i in range(n_reactions)}
            )
            lp.optimize()
            if lp.status != OPTIMAL:
                continue
            vector = np.array([var.primal for var in flux_vars], dtype=float)
            vector[np.abs(vector) < v_tol] = 0.0
            if not vector.any():
                continue
            candidates.append((int(np.count_nonzero(vector)), vector))

        if not candidates:
            # Neither direction admits anything independent of what is already
            # held, so the bounded null space is exhausted.
            break
        basis = np.column_stack([basis, min(candidates, key=lambda c: c[0])[1]])

    logger.debug(
        "Fast-SNP found %d of at most %d basis vectors", basis.shape[1], target
    )
    return basis

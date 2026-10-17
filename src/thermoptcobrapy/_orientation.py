"""Orienting reactions and taming bounds, as every ThermOpt routine does first.

Each entry point in the MATLAB package opens with the same two steps:

.. code-block:: matlab

    IrR = model.ub<=0;                  % irreversible in the reverse direction
    temp = model.ub(IrR);
    model.S(:,IrR) = -model.S(:,IrR);
    model.ub(IrR) = -1*model.lb(IrR);
    model.lb(IrR) = -1*temp;

    temp1 = max(abs(model.lb));         % normalise the extreme bounds to 1000
    temp2 = max(abs(model.ub));
    model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
    model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;

The first is a change of variables: negating a reaction's column and negating
its flux cancel out, so the network is unchanged and afterwards nothing is
restricted to negative flux. Every later step is then free to read "carries
flux" as "carries positive flux".

The second is narrower than it looks. It does not rescale the model: it finds
the single largest bound magnitude and clamps only the reactions that attain
it to 1000. Reactions at 999 are left alone. That is what the MATLAB does, so
that is what this does.
"""

from typing import TYPE_CHECKING, List, Tuple

if TYPE_CHECKING:
    from cobra.core import Model

#: The magnitude the extreme bounds are clamped to.
BOUND_NORM = 1000.0


def flip_reverse_only(model: "Model") -> List[str]:
    """Convert every reaction irreversible in reverse to forward, in place.

    Parameters
    ----------
    model : cobra.Model
        The model to reorient. Modified in place, so hand it a copy.

    Returns
    -------
    list of str
        The identifiers of the reactions that were flipped, in model order.

    Notes
    -----
    A reaction qualifies when its upper bound is at most zero and its lower
    bound is below zero, which is the MATLAB's ``model.ub<=0`` restricted to
    the reactions that can actually carry flux. A reaction pinned shut at
    ``(0, 0)`` is left alone: flipping it would be a no-op that still
    rewrites its stoichiometry.

    """
    flipped = []
    for reaction in model.reactions:
        if reaction.upper_bound > 0.0 or reaction.lower_bound >= 0.0:
            continue
        # Adding -2x the current coefficients leaves -x of each behind.
        reaction.add_metabolites(
            {met: -2.0 * coef for met, coef in reaction.metabolites.items()},
            combine=True,
        )
        lower, upper = reaction.bounds
        reaction.bounds = (-upper, -lower)
        flipped.append(reaction.id)
    return flipped


def normalise_bounds(model: "Model") -> None:
    """Clamp the largest bound magnitudes to ``BOUND_NORM``, in place.

    Parameters
    ----------
    model : cobra.Model
        The model to adjust. Modified in place.

    Notes
    -----
    Only the reactions attaining the maximum are touched, separately for the
    lower and the upper bounds, and a bound of exactly zero is never given a
    sign. Both bounds of a reaction are written together, so the model never
    passes through a state where the lower bound exceeds the upper one.

    """
    lowers = [rxn.lower_bound for rxn in model.reactions]
    uppers = [rxn.upper_bound for rxn in model.reactions]
    if not lowers:
        return
    max_lower = max(abs(value) for value in lowers)
    max_upper = max(abs(value) for value in uppers)

    for reaction in model.reactions:
        lower, upper = reaction.bounds
        if max_lower and abs(lower) == max_lower:
            lower = (1.0 if lower > 0 else -1.0) * BOUND_NORM
        if max_upper and abs(upper) == max_upper:
            upper = (1.0 if upper > 0 else -1.0) * BOUND_NORM
        if (lower, upper) != reaction.bounds:
            reaction.bounds = (min(lower, upper), max(lower, upper))


def prepare(model: "Model") -> Tuple["Model", List[str]]:
    """Return the working copy every ThermOpt entry point starts from.

    Parameters
    ----------
    model : cobra.Model
        The model to prepare. It is left unchanged.

    Returns
    -------
    tuple of (cobra.Model, list of str)
        A copy in which the reactions irreversible in the reverse direction
        have been converted to forward irreversible reactions and the extreme
        bounds clamped, plus the identifiers of the reactions that were
        converted. The copy is the ``modModel`` the MATLAB returns, and
        results refer to it rather than to the model that was passed in.

    """
    working = model.copy()
    flipped = flip_reverse_only(working)
    normalise_bounds(working)
    return working, flipped


def is_reversible(reaction) -> bool:
    """Return whether a reaction may carry negative flux.

    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to inspect, already oriented by :func:`flip_reverse_only`.

    Returns
    -------
    bool
        Whether the lower bound is negative, which is the MATLAB's
        ``model.lb < 0``.

    """
    return reaction.lower_bound < 0.0

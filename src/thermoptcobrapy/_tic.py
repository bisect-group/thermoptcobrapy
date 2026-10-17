"""Isolating the part of a model where thermodynamically infeasible cycles live.

Every ThermOpt routine narrows the problem the same way before doing anything
expensive. Close all the exchange reactions and drop any positive lower bound,
and the only flux the model can still carry is flux that goes round in a
circle: with nothing entering or leaving, a steady state has no other option.
The reactions that remain flux consistent under those conditions are exactly
the reactions that can take part in a cycle.

That is a consistency check, done with SPECTRA, mirroring the MATLAB's own
``spectraCC``.
"""

from typing import TYPE_CHECKING, List, Optional, Set

from ._subset import subset_model
from .exceptions import ThermOptError

if TYPE_CHECKING:
    from cobra.core import Model


def exchange_reaction_ids(model: "Model") -> Set[str]:
    """Return the identifiers of the model's exchange reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect.

    Returns
    -------
    set of str
        The reactions that touch exactly one metabolite, which is what
        MATLAB's ``findExcRxns`` detects and what cobrapy calls a boundary
        reaction. Exchanges, demands and sinks all qualify.

    """
    return {rxn.id for rxn in model.reactions if rxn.boundary}


def consistent_reaction_ids(
    model: "Model", tol: float, seed: Optional[int] = None
) -> Set[str]:
    """Return the flux consistent reactions of a model.

    Parameters
    ----------
    model : cobra.Model
        The model to check.
    tol : float
        The flux threshold a reaction has to reach.
    seed : int, optional
        A seed for SPECTRA's own randomisation (default None).

    Returns
    -------
    set of str
        The identifiers of the reactions that can carry flux.

    Raises
    ------
    ThermOptError
        If ``spectra_cobra`` is not installed. It is a required dependency,
        not an optional one: every consistency check in this package goes
        through it.

    """
    try:
        from spectra_cobra import consistent_reaction_ids as spectra_ids
    except ImportError as error:
        raise ThermOptError(
            "spectra_cobra is not installed, but every consistency check in "
            "this package depends on it. Clone and install it: "
            "https://github.com/bisect-group/spectra-cobra"
        ) from error

    ids, _ = spectra_ids(model, tol=tol, seed=seed)
    return ids


def cycle_submodel(model: "Model") -> "Model":
    """Return the model with its exchanges shut and positive lower bounds freed.

    Parameters
    ----------
    model : cobra.Model
        The model to restrict. It is left unchanged.

    Returns
    -------
    cobra.Model
        A copy in which every boundary reaction is pinned to zero flux and no
        reaction is forced to carry flux. Whatever stays consistent in it is
        cycling.

    """
    restricted = model.copy()
    for reaction in restricted.reactions:
        lower, upper = reaction.bounds
        if lower > 0.0:
            lower = 0.0
        if reaction.boundary:
            lower, upper = 0.0, 0.0
        reaction.bounds = (lower, upper)
    return restricted


def cycle_reaction_ids(
    model: "Model", tol: float, seed: Optional[int] = None
) -> List[str]:
    """Return the reactions that can take part in an internal cycle.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect, already oriented. It is left unchanged.
    tol : float
        The flux threshold a reaction has to reach to count as cycling.
    seed : int, optional
        A seed for SPECTRA's own randomisation (default None).

    Returns
    -------
    list of str
        The identifiers, in model order.

    """
    consistent = consistent_reaction_ids(cycle_submodel(model), tol, seed=seed)
    return [rxn.id for rxn in model.reactions if rxn.id in consistent]


def cycle_basis(model: "Model", tol: float, seed: Optional[int] = None):
    """Return the cycle reactions of a model and a sparse basis for them.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect, already oriented. It is left unchanged.
    tol : float
        The flux threshold a reaction has to reach to count as cycling.
    seed : int, optional
        A seed for the randomisation in both the consistency pre-check
        (SPECTRA's own) and Fast-SNP's random weights (default None). The
        same seed is used for both.

    Returns
    -------
    tuple of (list of str, numpy.ndarray)
        The cycling reactions in model order, and the ``n_tic x K`` basis
        whose rows follow them. Both are empty when the model has no cycles.

    """
    import numpy as np

    from ._nullspace import fast_snp

    restricted = cycle_submodel(model)
    consistent = consistent_reaction_ids(restricted, tol, seed=seed)
    tic_ids = [rxn.id for rxn in model.reactions if rxn.id in consistent]
    if not tic_ids:
        return [], np.zeros((0, 0))

    reduced = subset_model(restricted, consistent)
    # subset_model keeps the source order, so the basis rows line up with
    # tic_ids without any further bookkeeping.
    basis = fast_snp(reduced, v_tol=tol, seed=seed)
    assert [rxn.id for rxn in reduced.reactions] == tic_ids
    return tic_ids, basis

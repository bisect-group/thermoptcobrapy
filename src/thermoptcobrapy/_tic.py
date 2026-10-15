"""Isolating the part of a model where thermodynamically infeasible cycles live.

Every ThermOpt routine narrows the problem the same way before doing anything
expensive. Close all the exchange reactions and drop any positive lower bound,
and the only flux the model can still carry is flux that goes round in a
circle: with nothing entering or leaving, a steady state has no other option.
The reactions that remain flux consistent under those conditions are exactly
the reactions that can take part in a cycle.

That is a consistency check, which is why the MATLAB reaches for ``spectraCC``
and falls back to ``fastcc``. The same preference is kept here.
"""

from logging import getLogger
from typing import TYPE_CHECKING, List, Set

from ._subset import subset_model

if TYPE_CHECKING:
    from cobra.core import Model

logger = getLogger(__name__)


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


def consistent_reaction_ids(model: "Model", tol: float) -> Set[str]:
    """Return the flux consistent reactions of a model.

    Parameters
    ----------
    model : cobra.Model
        The model to check.
    tol : float
        The flux threshold a reaction has to reach.

    Returns
    -------
    set of str
        The identifiers of the reactions that can carry flux.

    Notes
    -----
    Prefers SPECTRA's consistency check when it is installed and falls back
    to cobrapy's ``fastcc``, mirroring the MATLAB's ``if exist('spectraCC')``.
    The two agree on which reactions are blocked; SPECTRA is the faster of
    them on larger models.

    """
    try:
        from spectra_cobra import consistent_reaction_ids as spectra_ids
    except ImportError:
        from cobra.flux_analysis import fastcc

        logger.debug("spectra_cobra is not installed; using fastcc")
        return {rxn.id for rxn in fastcc(model, flux_threshold=tol).reactions}

    ids, _ = spectra_ids(model, tol=tol)
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


def cycle_reaction_ids(model: "Model", tol: float) -> List[str]:
    """Return the reactions that can take part in an internal cycle.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect, already oriented. It is left unchanged.
    tol : float
        The flux threshold a reaction has to reach to count as cycling.

    Returns
    -------
    list of str
        The identifiers, in model order.

    """
    consistent = consistent_reaction_ids(cycle_submodel(model), tol)
    return [rxn.id for rxn in model.reactions if rxn.id in consistent]


def cycle_basis(model: "Model", tol: float, seed=None):
    """Return the cycle reactions of a model and a sparse basis for them.

    Parameters
    ----------
    model : cobra.Model
        The model to inspect, already oriented. It is left unchanged.
    tol : float
        The flux threshold a reaction has to reach to count as cycling.
    seed : int, optional
        A seed for Fast-SNP's random weights (default None).

    Returns
    -------
    tuple of (list of str, numpy.ndarray)
        The cycling reactions in model order, and the ``n_tic x K`` basis
        whose rows follow them. Both are empty when the model has no cycles.

    """
    import numpy as np

    from ._nullspace import fast_snp

    restricted = cycle_submodel(model)
    consistent = consistent_reaction_ids(restricted, tol)
    tic_ids = [rxn.id for rxn in model.reactions if rxn.id in consistent]
    if not tic_ids:
        return [], np.zeros((0, 0))

    reduced = subset_model(restricted, consistent)
    # subset_model keeps the source order, so the basis rows line up with
    # tic_ids without any further bookkeeping.
    basis = fast_snp(reduced, v_tol=tol, seed=seed)
    assert [rxn.id for rxn in reduced.reactions] == tic_ids
    return tic_ids, basis

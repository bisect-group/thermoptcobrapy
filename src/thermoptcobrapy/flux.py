"""Detecting and removing loops from a flux distribution, ThermOptFlux.

A flux distribution from FBA or from flux sampling satisfies the mass
balance, but may carry any amount of flux through a thermodynamically
infeasible cycle (TIC), because nothing in the problem penalises it. Given
the TICs of the model, this projects the distribution to the nearest one in
the thermodynamically feasible flux space.

A TIC is only subtracted when the distribution runs *every* one of its
reactions in the TIC's own direction. The amount removed is set by the
reaction carrying the least flux, so that one lands exactly on zero and no
reaction is pushed through it into the opposite direction.
"""

from logging import getLogger
from typing import Dict, Mapping, Sequence, Union

logger = getLogger(__name__)

FluxLike = Union[Mapping[str, float], "object"]


def _as_dict(fluxes: FluxLike) -> Dict[str, float]:
    """Return a plain dict of fluxes from a mapping or a pandas Series.

    Parameters
    ----------
    fluxes : mapping of {str: float} or pandas.Series
        The flux distribution.

    Returns
    -------
    dict of {str: float}
        The same values as a dict.

    """
    if hasattr(fluxes, "to_dict"):
        return dict(fluxes.to_dict())
    return dict(fluxes)


def thermopt_flux(
    fluxes: FluxLike,
    cycles: Sequence[Sequence[str]],
    directions: Sequence[Sequence[float]],
) -> Dict[str, float]:
    """Return the flux distribution with its TIC loops removed.

    Parameters
    ----------
    fluxes : mapping of {str: float} or pandas.Series
        The flux of every reaction, as FBA or sampling produced it.
    cycles : sequence of sequence of str
        The TICs, each a list of reaction identifiers. This is what
        :func:`~thermoptcobrapy.enumeration.thermopt_enum_milp` returns.
    directions : sequence of sequence of float
        The relative flux coefficient of each reaction in the matching TIC,
        element for element with `cycles`.

    Returns
    -------
    dict of {str: float}
        A new distribution. The input is not modified.

    Notes
    -----
    Reactions named in a TIC but absent from `fluxes` make that TIC
    inapplicable rather than an error, which is what happens when the TICs
    were enumerated on a GEM and the distribution comes from a sub-model.

    """
    corrected = _as_dict(fluxes)

    for cycle, direction in zip(cycles, directions):
        if len(cycle) != len(direction):
            raise ValueError(
                f"cycle of {len(cycle)} reactions carries {len(direction)} "
                "coefficients; they have to correspond element for element"
            )
        if any(rxn_id not in corrected for rxn_id in cycle):
            logger.debug("skipping a cycle whose reactions are not all present")
            continue

        present = [corrected[rxn_id] for rxn_id in cycle]
        running = sum(
            1
            for value, coefficient in zip(present, direction)
            if (value > 0 and coefficient > 0) or (value < 0 and coefficient < 0)
        )
        if running != len(cycle):
            continue

        # Scale the cycle so its smallest-carrying reaction is cancelled
        # exactly; subtracting more would reverse that reaction.
        smallest = min(range(len(cycle)), key=lambda i: abs(present[i]))
        scale = abs(present[smallest]) / abs(direction[smallest])
        for rxn_id, coefficient in zip(cycle, direction):
            corrected[rxn_id] -= coefficient * scale

    return corrected


def cycling_removed(before: FluxLike, after: Mapping[str, float]) -> Dict[str, float]:
    """Return how much flux each reaction lost to loop removal.

    Parameters
    ----------
    before : mapping of {str: float} or pandas.Series
        The distribution before correction.
    after : mapping of {str: float}
        The distribution after.

    Returns
    -------
    dict of {str: float}
        The difference per reaction, for the reactions that changed.

    """
    original = _as_dict(before)
    return {
        rxn_id: original[rxn_id] - after[rxn_id]
        for rxn_id in original
        if abs(original[rxn_id] - after.get(rxn_id, 0.0)) > 0.0
    }

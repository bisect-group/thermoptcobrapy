"""Splitting reversible reactions into a forward and a reverse half.

ThermOptCC needs this when it is asked for directions rather than a plain
blocked/unblocked verdict. Working on the split model lets each direction be
tested on its own: the search drops a reaction from the core as soon as it is
shown to carry flux, so without the split a reaction proved to run forward
would never be asked whether it can also run backward.

cobrapy used to ship ``convert_to_irreversible`` and no longer does, so the
part ThermOpt needs is kept here.
"""

from typing import TYPE_CHECKING, Dict, List, Tuple

if TYPE_CHECKING:
    from cobra.core import Model

#: Appended to a reaction's identifier to name its reverse half.
REVERSE_SUFFIX = "_ThermOpt_reverse"


def convert_to_irreversible(model: "Model") -> Tuple["Model", Dict[str, List[str]]]:
    """Return a copy with no reaction able to carry negative flux.

    Parameters
    ----------
    model : cobra.Model
        The model to split, already oriented so that nothing has a
        non-positive upper bound. It is left unchanged.

    Returns
    -------
    tuple of (cobra.Model, dict of {str: list of str})
        The split model, and a map from each original reaction to the halves
        that represent it: one identifier for a reaction that was already
        irreversible, two for one that was not, forward first.

    """
    from cobra.core import Reaction

    split = model.copy()
    mapping: Dict[str, List[str]] = {}
    new_halves = []

    for reaction in list(split.reactions):
        lower, upper = reaction.bounds
        if lower >= 0.0:
            mapping[reaction.id] = [reaction.id]
            continue

        reverse = Reaction(
            f"{reaction.id}{REVERSE_SUFFIX}",
            name=reaction.name,
            subsystem=reaction.subsystem,
            lower_bound=0.0,
            upper_bound=-lower,
        )
        reverse.add_metabolites(
            {met: -coef for met, coef in reaction.metabolites.items()}
        )
        reverse.gene_reaction_rule = reaction.gene_reaction_rule
        new_halves.append(reverse)

        reaction.bounds = (0.0, upper)
        mapping[reaction.id] = [reaction.id, reverse.id]

    if new_halves:
        split.add_reactions(new_halves)
    return split, mapping

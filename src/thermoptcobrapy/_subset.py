"""Building a sub-model from the reactions that are kept.

Deleting is quadratic: every reaction removed makes optlang re-index its
variable container and updates the solver once more. Building from what is
kept is linear and populates the solver in a single batch, which matters
because the context-specific routine discards most of a genome-scale model.
"""

from copy import deepcopy
from typing import TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    from cobra.core import Model


def subset_model(model: "Model", keep_ids: Iterable[str]) -> "Model":
    """Return a new model holding only the given reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to take the subset of. It is left unchanged.
    keep_ids : iterable of str
        The identifiers of the reactions to keep.

    Returns
    -------
    cobra.Model
        A new model carrying the kept reactions and everything attached to
        them, in the order the source had them.

    """
    from cobra import Model
    from cobra.core import Group
    from cobra.core.dictlist import DictList

    keep = set(keep_ids)
    reactions = [rxn for rxn in model.reactions if rxn.id in keep]
    metabolites = {met for rxn in reactions for met in rxn._metabolites}
    genes = {gene for rxn in reactions for gene in rxn._genes}

    # ``Reaction.__getstate__`` does not null ``_model``, unlike
    # ``Object.__getstate__``, so a plain deepcopy would drag the whole source
    # model and its solver along behind every reaction.
    detached = [(obj, obj._model) for obj in (*reactions, *metabolites, *genes)]
    for obj, _ in detached:
        obj._model = None
    try:
        copies = deepcopy(reactions)
    finally:
        for obj, owner in detached:
            obj._model = owner

    subset = Model(model.id)
    subset.name = model.name
    subset.compartments = dict(model.compartments)
    subset.annotation = deepcopy(model.annotation)
    subset.notes = deepcopy(model.notes)
    subset.add_reactions(copies)

    # The GPR is serialised as a string, so the genes are rebuilt by parsing
    # it and arrive carrying nothing but an identifier.
    for gene in subset.genes:
        original = model.genes.get_by_id(gene.id)
        gene.name = original.name
        gene.annotation = deepcopy(original.annotation)
        gene.notes = deepcopy(original.notes)

    for attribute in ("metabolites", "genes"):
        position = {obj.id: i for i, obj in enumerate(getattr(model, attribute))}
        current = list(getattr(subset, attribute))
        setattr(
            subset,
            attribute,
            DictList(sorted(current, key=lambda obj: position[obj.id])),
        )

    groups = []
    for group in model.groups:
        members = []
        for member in group.members:
            for container in (subset.reactions, subset.metabolites, subset.genes):
                if member.id in container:
                    members.append(container.get_by_id(member.id))
                    break
        new_group = Group(group.id, name=group.name, members=members, kind=group.kind)
        new_group.annotation = deepcopy(group.annotation)
        new_group.notes = deepcopy(group.notes)
        groups.append(new_group)
    if groups:
        subset.add_groups(groups)

    subset.solver = model.solver.interface
    subset.tolerance = model.tolerance
    coefficients = {
        subset.reactions.get_by_id(rxn.id): rxn.objective_coefficient
        for rxn in reactions
        if rxn.objective_coefficient
    }
    if coefficients:
        subset.objective = coefficients
    subset.objective.direction = model.objective.direction
    return subset

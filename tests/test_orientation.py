"""Test the orientation and bound-normalising step."""

from thermoptcobrapy import flip_reverse_only, normalise_bounds, prepare


def test_flip_turns_a_reverse_only_reaction_forward(reverse_only_model):
    """Test that a reaction with ub <= 0 is negated and its bounds swapped."""
    before = {
        met.id: coef
        for met, coef in reverse_only_model.reactions.R1.metabolites.items()
    }
    flipped = flip_reverse_only(reverse_only_model)

    assert flipped == ["R1"]
    assert reverse_only_model.reactions.R1.bounds == (0.0, 1000.0)
    after = {
        met.id: coef
        for met, coef in reverse_only_model.reactions.R1.metabolites.items()
    }
    assert after == {key: -value for key, value in before.items()}


def test_flip_leaves_a_forward_reaction_alone(cycle_model):
    """Test that nothing is flipped when no reaction is reverse-only."""
    before = {rxn.id: rxn.bounds for rxn in cycle_model.reactions}
    assert flip_reverse_only(cycle_model) == []
    assert {rxn.id: rxn.bounds for rxn in cycle_model.reactions} == before


def test_flip_skips_a_reaction_pinned_shut(cycle_model):
    """Test that a reaction fixed at zero is not rewritten.

    Its upper bound is zero, so a literal reading of ``ub <= 0`` would flip
    it, negating its stoichiometry to no purpose.

    """
    cycle_model.reactions.R4.bounds = (0.0, 0.0)
    before = {
        met.id: coef for met, coef in cycle_model.reactions.R4.metabolites.items()
    }
    assert "R4" not in flip_reverse_only(cycle_model)
    assert {
        met.id: coef for met, coef in cycle_model.reactions.R4.metabolites.items()
    } == before


def test_normalise_clamps_only_the_extremes(cycle_model):
    """Test that only the reactions attaining the maximum are clamped."""
    cycle_model.reactions.R1.bounds = (-5000.0, 5000.0)
    cycle_model.reactions.R2.bounds = (-10.0, 10.0)
    normalise_bounds(cycle_model)

    assert cycle_model.reactions.R1.bounds == (-1000.0, 1000.0)
    # R2 is nowhere near the maximum and must be untouched.
    assert cycle_model.reactions.R2.bounds == (-10.0, 10.0)


def test_prepare_leaves_the_input_alone(reverse_only_model):
    """Test that preparing returns a copy rather than editing in place."""
    before = reverse_only_model.reactions.R1.bounds
    working, flipped = prepare(reverse_only_model)

    assert flipped == ["R1"]
    assert reverse_only_model.reactions.R1.bounds == before
    assert working.reactions.R1.bounds != before

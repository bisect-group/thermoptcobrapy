"""Test the removal of cycling from a flux distribution."""

import pytest

from thermoptcobrapy import cycling_removed, thermopt_flux


def test_a_fully_cycling_distribution_is_zeroed():
    """Test that flux running entirely round a cycle is removed."""
    fluxes = {"R1": 5.0, "R2": 5.0, "R3": 5.0, "R4": 2.0}
    corrected = thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])

    assert corrected == {"R1": 0.0, "R2": 0.0, "R3": 0.0, "R4": 2.0}


def test_only_the_cycling_part_is_removed():
    """Test that the correction stops at the least-carrying reaction.

    Taking more would push that reaction through zero into the opposite
    direction, which would be inventing flux rather than removing it.

    """
    fluxes = {"R1": 7.0, "R2": 5.0, "R3": 6.0}
    corrected = thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])

    assert corrected == {"R1": 2.0, "R2": 0.0, "R3": 1.0}


def test_a_cycle_not_being_run_is_left_alone():
    """Test that a cycle only counts when every reaction follows it."""
    fluxes = {"R1": 5.0, "R2": 5.0, "R3": -5.0}
    corrected = thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])

    assert corrected == fluxes


def test_direction_signs_are_respected():
    """Test that a cycle running a reaction backwards is handled."""
    fluxes = {"R1": 4.0, "R2": -4.0}
    corrected = thermopt_flux(fluxes, [["R1", "R2"]], [[1.0, -1.0]])

    assert corrected == {"R1": 0.0, "R2": 0.0}


def test_the_input_is_not_modified():
    """Test that the distribution handed in is left as it was."""
    fluxes = {"R1": 5.0, "R2": 5.0, "R3": 5.0}
    thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])

    assert fluxes == {"R1": 5.0, "R2": 5.0, "R3": 5.0}


def test_a_cycle_naming_absent_reactions_is_skipped():
    """Test that a cycle from another model does not raise."""
    fluxes = {"R1": 5.0}
    assert thermopt_flux(fluxes, [["R1", "R9"]], [[1.0, 1.0]]) == {"R1": 5.0}


def test_mismatched_lengths_are_rejected():
    """Test that a cycle and its coefficients have to correspond."""
    with pytest.raises(ValueError, match="element for element"):
        thermopt_flux({"R1": 1.0}, [["R1"]], [[1.0, 1.0]])


def test_cycling_removed_reports_the_difference():
    """Test the helper that says how much each reaction lost."""
    fluxes = {"R1": 5.0, "R2": 5.0, "R3": 5.0, "R4": 2.0}
    corrected = thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])
    removed = cycling_removed(fluxes, corrected)

    assert removed == {"R1": 5.0, "R2": 5.0, "R3": 5.0}

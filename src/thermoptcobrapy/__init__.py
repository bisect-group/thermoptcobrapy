"""ThermOpt for cobrapy.

Flux balance analysis constrains a network with the mass balance alone, and
the mass balance does not account for thermodynamically infeasible cycles
(TICs). Analogous to perpetual motion machines, these cycles violate the
second law of thermodynamics by cycling metabolites indefinitely without any
real change.

TICs lead genome-scale metabolic models (GEMs) to predict phenotypes that are
thermodynamically infeasible, and such phenotypes carry no biological
interpretation. FBA, FVA and flux sampling alike often predict maximum flux
through reactions involved in TICs, and TICs are inherited by every
context-specific model built on a network that still contains them.

This package offers four ways at them:

1. :func:`thermopt_cc` determines the thermodynamically feasible flux
   directions of every reaction, and thereby the stoichiometrically and
   thermodynamically blocked reactions.
2. :func:`thermopt_enum_milp` enumerates all TICs of a model and the flux
   direction of the reactions involved.
3. :func:`thermopt_flux` detects and removes loops from a flux distribution,
   projecting it to the nearest distribution in the thermodynamically
   feasible flux space.
4. :func:`thermoptics` constructs a thermodynamically consistent
   context-specific model (CSM).

The mixed-integer routines need a real mixed-integer solver. Gurobi and CPLEX
both qualify; GLPK does not, and is rejected rather than trusted when it
returns a solution its own constraints forbid.

Reference
---------
Pavan Kumar S. and Nirav Pravinbhai Bhatt (2025). ThermOptCobra:
Thermodynamically optimal construction and analysis of metabolic networks for
reliable phenotype predictions. iScience 28(8), 113005.
https://doi.org/10.1016/j.isci.2025.113005
"""

from ._irreversible import convert_to_irreversible
from ._nullspace import fast_snp
from ._orientation import flip_reverse_only, normalise_bounds, prepare
from ._subset import subset_model
from ._tic import cycle_basis, cycle_reaction_ids, cycle_submodel
from .consistency import (
    BLOCKED,
    FORWARD,
    REVERSE,
    REVERSIBLE,
    UNBLOCKED,
    thermodynamically_blocked_ids,
    thermopt_cc,
)
from .context import thermoptics
from .enumeration import reaction_order, thermopt_enum_milp
from .exceptions import (
    ThermOptError,
    ThermOptInfeasibleCoreError,
    ThermOptSolverError,
    ThermOptTimeLimitError,
)
from .flux import cycling_removed, thermopt_flux

__version__ = "0.1.0.dev0"

__all__ = (
    # Consistency.
    "thermopt_cc",
    "thermodynamically_blocked_ids",
    "BLOCKED",
    "UNBLOCKED",
    "FORWARD",
    "REVERSE",
    "REVERSIBLE",
    # Cycle enumeration.
    "thermopt_enum_milp",
    "reaction_order",
    # Flux correction.
    "thermopt_flux",
    "cycling_removed",
    # Context-specific extraction.
    "thermoptics",
    # Shared machinery, useful on its own.
    "fast_snp",
    "cycle_basis",
    "cycle_reaction_ids",
    "cycle_submodel",
    "prepare",
    "flip_reverse_only",
    "normalise_bounds",
    "convert_to_irreversible",
    "subset_model",
    # Errors.
    "ThermOptError",
    "ThermOptSolverError",
    "ThermOptInfeasibleCoreError",
    "ThermOptTimeLimitError",
)

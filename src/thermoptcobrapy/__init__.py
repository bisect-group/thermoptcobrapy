"""ThermOpt for cobrapy.

Flux balance analysis constrains a network with the mass balance alone, and
the mass balance has no opinion about flux that goes round in a circle. A set
of reactions whose net effect is nothing can carry arbitrary flux at steady
state. Thermodynamics forbids it: the free energy would have to fall all the
way round the cycle and arrive back where it started.

These *thermodynamically infeasible cycles* are not a curiosity. They inflate
maximum fluxes, make blocked reactions look usable, and are inherited by every
context-specific model built on top of a network that still contains them.

This package offers four ways at them:

1. :func:`thermopt_cc` reports, for every reaction, which directions it can
   carry flux in once cycling is excluded.
2. :func:`thermopt_enum_milp` lists the cycles themselves, with the relative
   flux through each reaction.
3. :func:`thermopt_flux` takes the cycling back out of a flux distribution
   that already has it.
4. :func:`thermoptics` extracts a context-specific model that is
   thermodynamically feasible by construction.

The mixed-integer routines need a real mixed-integer solver. Gurobi and CPLEX
both qualify; GLPK does not, and is rejected rather than trusted when it
returns a solution its own constraints forbid.

Reference
---------
Thermodynamically constrained reconstruction and analysis.
https://www.sciencedirect.com/science/article/pii/S2589004225012660
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

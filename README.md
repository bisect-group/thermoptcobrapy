# thermoptcobrapy

**ThermOpt for cobrapy.**

A Python port of the [ThermOptCOBRA](https://github.com/bisect-group/Protocol_ThermOptCOBRA)
MATLAB package, built on [cobrapy](https://github.com/opencobra/cobrapy).

Flux balance analysis constrains a network with the mass balance alone, and
the mass balance does not account for *thermodynamically infeasible cycles*
(TICs). Analogous to perpetual motion machines, these cycles violate the
second law of thermodynamics by cycling metabolites indefinitely without any
real change.

TICs lead genome-scale metabolic models (GEMs) to predict phenotypes that are
thermodynamically infeasible, and such phenotypes carry no biological
interpretation. FBA, FVA and flux sampling alike often predict maximum flux
through reactions involved in TICs, and TICs are inherited by every
context-specific model built on a network that still contains them.

ThermOpt addresses them in four ways:

| Routine | What it does |
|---|---|
| `thermopt_cc` | Determines the thermodynamically feasible flux directions of every reaction, and thereby the stoichiometrically and thermodynamically blocked reactions |
| `thermopt_enum_milp` | Enumerates all TICs of a model and the flux direction of the reactions involved |
| `thermopt_flux` | Detects and removes loops from a flux distribution, projecting it to the nearest distribution in the thermodynamically feasible flux space |
| `thermoptics` | Constructs a thermodynamically consistent context-specific model (CSM) |

Reference: Pavan Kumar S. and Nirav Pravinbhai Bhatt (2025). ThermOptCobra:
Thermodynamically optimal construction and analysis of metabolic networks for
reliable phenotype predictions. *iScience* **28**(8), 113005.
<https://doi.org/10.1016/j.isci.2025.113005>

## Installation

Every internal consistency check this package runs goes through
[spectra-cobra](https://github.com/bisect-group/spectra-cobra), so it is a
required dependency rather than an optional one. It is not published on
PyPI, so clone and install it first:

```bash
git clone https://github.com/bisect-group/spectra-cobra.git
cd spectra-cobra
pip install -e .
cd ..
```

Then thermoptcobrapy itself:

```bash
git clone https://github.com/bisect-group/thermoptcobrapy.git
cd thermoptcobrapy
pip install -e .
```

The mixed-integer routines need Gurobi or CPLEX:

```bash
pip install -e ".[gurobi]"
```

GLPK is not merely slow on these problems, it is wrong: it has been seen
reporting `optimal` with a reaction's indicator set while that reaction
carries no flux. Every solution is checked against its own indicator
constraints and a solver that returns one it forbids raises rather than
being trusted.


## Choosing `tol`

`tol` is the flux a reaction must reach to count as carrying flux, and it is
not free to be as small as you like. The indicator constraints are enforced
through terms like `ub_i * z_i`, so a binary sitting `integrality` away from
zero buys `ub_i * integrality` of flux for nothing. Once that reaches `tol`
the indicators stop meaning anything.

Gurobi will not accept an integrality tolerance below `1e-9`, so on a model
whose bounds reach 1000 the smallest workable threshold is about `5e-5`. The
package checks this up front and tells you what to use rather than returning
a quietly wrong answer.

Two consequences worth knowing:

- Reactions that cannot carry `tol` are reported blocked. That is the correct
  result at that threshold, not a failure.
- A model with a **forced** lower bound — a biomass reaction pinned to a
  growth rate, say — can be infeasible rather than merely restrictive. Biomass
  at 25.03 on an *E. coli* model demands trace-cofactor fluxes around `6e-5`
  (biotin's biomass coefficient is `2e-6`). That still works, but only in a
  narrow window: `tol = 5e-5` resolves it, `tol = 1e-4` is too coarse and
  makes the whole program infeasible. Releasing such bounds before asking
  which reactions *can* carry flux avoids the squeeze, and is anyway a
  different question from what happens at a fixed growth rate.

A note on `feasTol`: the MATLAB sets the solver's feasibility tolerance to
`1e-9`. Doing that with Gurobi 13 makes it return provably wrong optima — on
one model it proved "no reaction can be active" while a solution with 561
active reactions is feasible to `1e-12`, verified outside the solver. The
feasibility tolerance is therefore held at or above `1e-6` for these solves,
and your model's own tolerance is restored afterwards.

TIC enumeration
===============

Enumerates all thermodynamically infeasible cycles (TICs) of a model, with
the relative flux coefficient of each reaction involved.

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_enum_milp

   cycles, directions, cycle_rxns, oriented = thermopt_enum_milp(model)

   for cycle, direction in zip(cycles, directions):
       print(cycle, direction)

``cycles[i]`` is a list of reaction identifiers and ``directions[i]`` gives
the relative flux coefficient of each of them, element for element. Those
two are exactly what :doc:`thermopt_flux` consumes.

Where the TICs are looked for
-----------------------------

The enumerator begins by closing the model: every exchange reaction is pinned
shut and every positive lower bound released. With no input or output of
nutrients, any non-zero flux a steady state can carry is cycling flux, so
whatever remains flux consistent under those conditions is exactly what can
take part in a TIC. That is a flux consistency check, done with SPECTRA.

The model handed back has its exchanges *open* — it is the oriented model
the TICs refer to (the MATLAB's ``modModel``), not the closed one they were
found in.

How it works
------------

``thermopt_enum_milp`` minimises the *number* of active reactions, which
yields the TIC itself rather than a TIC plus whatever else the solver left
switched on. It forces one reaction to carry flux, solves for the smallest
TIC through it, and excludes that TIC with a no-good cut so the next solve
must produce a different one. Each TIC found also narrows the model the next
solve runs on to that TIC's own reactions, which is what keeps the repeated
solves cheap.

It runs at a flux threshold of 1. In a closed network the scale of a TIC is
arbitrary; only the set of reactions involved and their relative flux
coefficients matter. It needs a real mixed-integer solver: Gurobi or CPLEX.

A worked example
----------------

A network of five reversible internal reactions over four metabolites, plus
two exchange reactions. ``R1``, ``R2`` and ``R3`` form the cycle
``A → B → C → A``, and ``R5`` closes a second one directly between ``B`` and
``A``:

.. list-table::
   :header-rows: 1
   :widths: 16 34 50

   * - Reaction
     - Stoichiometry
     - Bounds
   * - ``R1``
     - ``A → B``
     - (-1000, 1000)
   * - ``R2``
     - ``B → C``
     - (-1000, 1000)
   * - ``R3``
     - ``C → A``
     - (-1000, 1000)
   * - ``R4``
     - ``A → D``
     - (-1000, 1000)
   * - ``R5``
     - ``B → A``
     - (-1000, 1000)
   * - ``EX_A``
     - ``A ↔``
     - (-1000, 1000)
   * - ``EX_D``
     - ``D ↔``
     - (-1000, 1000)

``R4`` and the two exchange reactions take part in no TIC: ``R4`` drains
``A`` to ``D``, which only ``EX_D`` touches, so none of the three can carry
flux once the exchanges are shut. :doc:`../tutorials/cycles` builds this
network in cobrapy if you want to run it.

Three TICs exist, each feasible in two directions, and the enumerator finds
all six:

.. code-block:: text

   ['R1', 'R5']              [ 1.0,  1.0]
   ['R1', 'R2', 'R3']        [ 1.0,  1.0,  1.0]
   ['R1', 'R5']              [-1.0, -1.0]
   ['R1', 'R2', 'R3']        [-1.0, -1.0, -1.0]
   ['R2', 'R3', 'R5']        [ 1.0,  1.0, -1.0]
   ['R2', 'R3', 'R5']        [-1.0, -1.0,  1.0]

Reading the first row: ``R1`` (``A → B``) and ``R5`` (``B → A``) both run
forward, carrying equal flux, which returns ``A`` to itself. The second is
the three-reaction cycle ``A → B → C → A``.

The third pair is the composite, and shows what the sign of a coefficient
means: ``R2`` and ``R3`` run forward while ``R5`` runs in reverse
(``A → B``), so ``A → B`` via ``R5`` reversed, then ``B → C → A`` via ``R2``
and ``R3``. A coefficient of ``-1.0`` is the reaction running against the
direction its stoichiometry is written in, not a different reaction.

Note that a TIC made entirely of reversible reactions is also feasible in
the reverse direction, so its mirror image is recorded rather than searched
for — which is why each of the three appears twice.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_enum_milp

.. autofunction:: thermoptcobrapy.reaction_order

See also
--------

* :doc:`../tutorials/cycles` — worked examples
* :doc:`thermopt_flux` — what to do with the TICs once you have them

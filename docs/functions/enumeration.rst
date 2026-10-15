Cycle enumeration
=================

Lists the thermodynamically infeasible cycles of a model, with the relative
flux through each of their reactions.

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_enum_lp

   cycles, directions, cycle_rxns, oriented = thermopt_enum_lp(model)

   for cycle, direction in zip(cycles, directions):
       print(cycle, direction)

``cycles[i]`` is a list of reaction identifiers and ``directions[i]`` gives
the flux through each of them, element for element. Those two are exactly
what :doc:`thermopt_flux` consumes.

Where the cycles are looked for
-------------------------------

Both enumerators begin by closing the model: every exchange reaction is
pinned shut and every positive lower bound released. With nothing entering or
leaving, the only flux a steady state can carry is flux that goes round in a
circle, so whatever stays flux consistent under those conditions is exactly
what can take part in a cycle. That is a consistency check, which is why the
implementation reaches for SPECTRA and falls back to ``fastcc``.

The model handed back has its exchanges *open* — it is the oriented model the
cycles refer to, not the closed one they were found in.

Two enumerators
---------------

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * -
     - ``thermopt_enum_lp``
     - ``thermopt_enum_milp``
   * - Objective
     - smallest one-norm of flux
     - fewest active reactions
   * - Cost
     - lower
     - higher
   * - Cycles returned
     - may carry passengers
     - tighter

``thermopt_enum_lp`` forces one reaction to carry flux and minimises
:math:`\|v\|_1`, which lands on a short cycle through it. It then blocks each
of that cycle's other reactions in turn and repeats, so cycles sharing
reactions are all reached; the recursion is memoised on the set of blocked
reactions.

``thermopt_enum_milp`` minimises the *number* of active reactions, which is
the cycle one actually wants, and excludes each cycle it finds with a no-good
cut so the next solve must produce a different one.

Both run at a flux threshold of 1. In a closed network the scale of a cycle
is arbitrary; only its support matters.

A worked example
----------------

Take ``A→B→C→A`` with an extra ``B→A``. There are three cycles, each in two
directions, and both enumerators find all six:

.. code-block:: text

   ['R1', 'R5']              [ 1.0,  1.0]
   ['R1', 'R2', 'R3']        [ 1.0,  1.0,  1.0]
   ['R1', 'R5']              [-1.0, -1.0]
   ['R1', 'R2', 'R3']        [-1.0, -1.0, -1.0]
   ['R2', 'R3', 'R5']        [ 1.0,  1.0, -1.0]
   ['R2', 'R3', 'R5']        [-1.0, -1.0,  1.0]

The third pair is the composite: ``R2`` and ``R3`` forward, against ``R5``
reversed. Note that a cycle made entirely of reversible reactions runs
backwards too, so its mirror image is recorded rather than searched for.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_enum_lp

.. autofunction:: thermoptcobrapy.thermopt_enum_milp

.. autofunction:: thermoptcobrapy.reaction_order

See also
--------

* :doc:`../tutorials/cycles` — worked examples
* :doc:`thermopt_flux` — what to do with the cycles once you have them

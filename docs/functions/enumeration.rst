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

Take ``A→B→C→A`` with an extra ``B→A``. There are three TICs, each feasible
in two directions, and the enumerator finds all six:

.. code-block:: text

   ['R1', 'R5']              [ 1.0,  1.0]
   ['R1', 'R2', 'R3']        [ 1.0,  1.0,  1.0]
   ['R1', 'R5']              [-1.0, -1.0]
   ['R1', 'R2', 'R3']        [-1.0, -1.0, -1.0]
   ['R2', 'R3', 'R5']        [ 1.0,  1.0, -1.0]
   ['R2', 'R3', 'R5']        [-1.0, -1.0,  1.0]

The third pair is the composite: ``R2`` and ``R3`` forward, against ``R5``
reversed. Note that a TIC made entirely of reversible reactions is also
feasible in the reverse direction, so its mirror image is recorded rather
than searched for.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_enum_milp

.. autofunction:: thermoptcobrapy.reaction_order

See also
--------

* :doc:`../tutorials/cycles` — worked examples
* :doc:`thermopt_flux` — what to do with the TICs once you have them

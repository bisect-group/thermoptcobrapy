thermopt_flux
=============

Detects and removes loops from a flux distribution, projecting it to the
nearest distribution in the thermodynamically feasible flux space.

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_enum_milp, thermopt_flux

   cycles, directions, _, oriented = thermopt_enum_milp(model)
   solution = oriented.optimize()

   corrected = thermopt_flux(solution.fluxes, cycles, directions)

The input is left unmodified; a new dict is returned. A pandas ``Series``,
as ``solution.fluxes`` is, works directly.

How it works
------------

Take the network from :doc:`../tutorials/cycles` — four metabolites, four
internal reactions and two exchange reactions. All six are reversible,
bounded ``(-1000, 1000)``, so the arrow below is only the reference
direction, the one a *positive* flux runs in:

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Reaction
     - Stoichiometry
     - Role
   * - ``R1``
     - ``A ↔ B``
     - in the TIC
   * - ``R2``
     - ``B ↔ C``
     - in the TIC
   * - ``R3``
     - ``C ↔ A``
     - in the TIC
   * - ``R4``
     - ``A ↔ D``
     - throughput
   * - ``EX_A``
     - ``A ↔``
     - throughput, takes ``A`` up
   * - ``EX_D``
     - ``D ↔``
     - throughput, exports ``D``

``R1``, ``R2`` and ``R3`` are the TIC: together they return ``A`` to itself
and produce nothing. The other three carry genuine throughput — ``A`` enters
through ``EX_A`` and leaves as ``D`` through ``EX_D``.

Now an FBA solution on that network. It converts 10 units of ``A`` into
``D``, and carries 3 units round the TIC on top:

.. code-block:: python

   >>> fluxes = {"R1": 3.0, "R2": 3.0, "R3": 3.0,
   ...           "R4": 10.0, "EX_A": -10.0, "EX_D": 10.0}

Every metabolite balances, so FBA has no objection: ``B`` and ``C`` are
produced and consumed by the TIC alone, and the 3 units on ``R1``, ``R2``
and ``R3`` cancel out at ``A``. The 3 units are nonetheless going round in a
circle and producing nothing.

Hand the TIC and its relative flux coefficients to ``thermopt_flux`` and it
goes:

.. code-block:: python

   >>> corrected = thermopt_flux(fluxes, [["R1", "R2", "R3"]],
   ...                            [[1.0, 1.0, 1.0]])
   >>> corrected
   {'R1': 0.0, 'R2': 0.0, 'R3': 0.0,
    'R4': 10.0, 'EX_A': -10.0, 'EX_D': 10.0}

   >>> cycling_removed(fluxes, corrected)
   {'R1': 3.0, 'R2': 3.0, 'R3': 3.0}

The throughput is untouched and every metabolite still balances. Only the
loop is gone.

Two rules decide what gets subtracted
-------------------------------------

**A TIC is subtracted only when every one of its reactions runs in the
TIC's own direction.** Add ``R5`` (``B → A``), which closes a second TIC
with ``R1``, and consider a solution that loads only that one:

.. code-block:: python

   >>> fluxes = {"R1": 4.0, "R2": 0.0, "R3": 0.0, "R5": 4.0,
   ...           "R4": 10.0, "EX_A": -10.0, "EX_D": 10.0}
   >>> thermopt_flux(fluxes, [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])
   {'R1': 4.0, 'R2': 0.0, 'R3': 0.0, 'R5': 4.0,
    'R4': 10.0, 'EX_A': -10.0, 'EX_D': 10.0}

Nothing changed: ``R2`` and ``R3`` carry no flux, so the ``R1``-``R2``-``R3``
TIC is not present in this distribution and subtracting it would invent flux
rather than remove it. The ``R1``-``R5`` TIC *is* present, and is removed:

.. code-block:: python

   >>> thermopt_flux(fluxes, [["R1", "R5"]], [[1.0, 1.0]])
   {'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R5': 0.0,
    'R4': 10.0, 'EX_A': -10.0, 'EX_D': 10.0}

**The amount removed is set by the reaction carrying the least flux**, so
that one lands exactly on zero and no reaction is pushed through zero into
the opposite direction. With both TICs loaded unequally — 4 units round
``R1``-``R5`` and 3 more round ``R1``-``R2``-``R3``, so ``R1`` carries 7:

.. code-block:: python

   >>> fluxes = {"R1": 7.0, "R2": 3.0, "R3": 3.0, "R5": 4.0,
   ...           "R4": 10.0, "EX_A": -10.0, "EX_D": 10.0}
   >>> thermopt_flux(fluxes, [["R1", "R5"]], [[1.0, 1.0]])
   {'R1': 3.0, 'R2': 3.0, 'R3': 3.0, 'R5': 0.0,
    'R4': 10.0, 'EX_A': -10.0, 'EX_D': 10.0}

``R5`` carried the least of the two at 4, so 4 is what came off both. ``R1``
keeps the 3 units belonging to the other TIC, which the next TIC in the list
then removes. Passing all six TICs the enumerator returns for this network
clears every loop in one call and still leaves ``R4``, ``EX_A`` and ``EX_D``
at 10, -10 and 10.

On a genome-scale model
-----------------------

Most FBA optima carry no loop, and on those this is a no-op — which is the
correct result, not a failure. Drive one deliberately and the effect is
unmistakable. Maximising ``ADK1`` (adenylate kinase, a reaction well known
to participate in a TIC) on a 1523-reaction *E. coli* model:

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * -
     - Before
     - After
   * - ``ADK1`` flux
     - 1000
     - 0
   * - Total :math:`\|v\|_1`
     - 3000
     - 0
   * - :math:`\max |S v|`
     - 0
     - 0

The whole distribution was loop flux, and the correction removes it while
the mass balance continues to hold.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_flux

.. autofunction:: thermoptcobrapy.cycling_removed

See also
--------

* :doc:`enumeration` — where ``cycles`` and ``directions`` come from
* :doc:`thermopt_cc` — the blocked reactions behind the TICs

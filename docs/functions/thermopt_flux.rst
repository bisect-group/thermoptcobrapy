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

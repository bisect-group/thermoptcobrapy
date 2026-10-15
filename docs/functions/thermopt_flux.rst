thermopt_flux
=============

Takes the infeasible cycling back out of a flux distribution that already
has it.

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_enum_lp, thermopt_flux

   cycles, directions, _, oriented = thermopt_enum_lp(model)
   solution = oriented.optimize()

   corrected = thermopt_flux(solution.fluxes, cycles, directions)

The input is left unmodified; a new dict is returned. A pandas ``Series``,
as ``solution.fluxes`` is, works directly.

How it works
------------

A flux vector from FBA or from sampling satisfies the mass balance but may
carry any amount of flux round a cycle, because nothing in the problem
charges for it. For each cycle, if the distribution runs **every** one of its
reactions the way the cycle does, a scaled copy of the cycle is subtracted.

The amount removed is set by the reaction carrying the least flux, so that
one lands exactly on zero. Taking more would push it through zero into the
opposite direction, which would be inventing flux rather than removing it:

.. code-block:: python

   >>> thermopt_flux({"R1": 7.0, "R2": 5.0, "R3": 6.0},
   ...               [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])
   {'R1': 2.0, 'R2': 0.0, 'R3': 1.0}

A cycle whose reactions are not all running its way is left alone entirely:

.. code-block:: python

   >>> thermopt_flux({"R1": 5.0, "R2": 5.0, "R3": -5.0},
   ...               [["R1", "R2", "R3"]], [[1.0, 1.0, 1.0]])
   {'R1': 5.0, 'R2': 5.0, 'R3': -5.0}

What it looks like on a real model
----------------------------------

Most FBA optima do not load a cycle, and on those this is a no-op — which is
the correct answer, not a failure. Drive one deliberately and the effect is
unmistakable. Maximising ``ADK1`` (adenylate kinase, a textbook infeasible
cycle) on a 1523-reaction *E. coli* model:

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

The distribution was pure cycling, and the correction removes it without
disturbing the mass balance.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_flux

.. autofunction:: thermoptcobrapy.cycling_removed

See also
--------

* :doc:`enumeration` — where ``cycles`` and ``directions`` come from

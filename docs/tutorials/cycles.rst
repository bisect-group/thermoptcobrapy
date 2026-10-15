Finding the cycles
==================

A network small enough to check by hand
---------------------------------------

Three reactions in a ring, plus a reaction that drains the ring to somewhere
useful:

.. code-block:: python

   from cobra import Model, Metabolite, Reaction

   model = Model("ring")
   A, B, C, D = (Metabolite(x) for x in "ABCD")

   reactions = []
   for rid, stoichiometry in (("R1", {A: -1, B: 1}),
                              ("R2", {B: -1, C: 1}),
                              ("R3", {C: -1, A: 1}),
                              ("R4", {A: -1, D: 1})):
       r = Reaction(rid, lower_bound=-1000, upper_bound=1000)
       r.add_metabolites(stoichiometry)
       reactions.append(r)
   for rid, met in (("EX_A", A), ("EX_D", D)):
       r = Reaction(rid, lower_bound=-1000, upper_bound=1000)
       r.add_metabolites({met: -1})
       reactions.append(r)
   model.add_reactions(reactions)
   model.solver = "gurobi"

``B`` and ``C`` are touched by nothing but the ring, so at steady state
``R1``, ``R2`` and ``R3`` must all carry the same flux — and that flux is a
cycle. FVA does not mind:

.. code-block:: python

   >>> from cobra.flux_analysis import flux_variability_analysis
   >>> flux_variability_analysis(model, fraction_of_optimum=0.0,
   ...                           processes=1).loc["R1"]
   minimum   -1000.0
   maximum    1000.0

Enumerating it
--------------

.. code-block:: python

   >>> from thermoptcobrapy import thermopt_enum_lp
   >>> cycles, directions, cycle_rxns, oriented = thermopt_enum_lp(model)
   >>> cycle_rxns
   ['R1', 'R2', 'R3']
   >>> for c, d in zip(cycles, directions):
   ...     print(c, d)
   ['R1', 'R2', 'R3'] [1.0, 1.0, 1.0]
   ['R1', 'R2', 'R3'] [-1.0, -1.0, -1.0]

Both directions are reported, because both are cycles. ``R4`` and the
exchanges are absent: they are not part of any ring.

Checking the answer
-------------------

A reported cycle should be a genuine steady-state mode of the closed network.
That is worth verifying rather than trusting:

.. code-block:: python

   import numpy as np
   from cobra.util.array import create_stoichiometric_matrix
   from thermoptcobrapy.enumeration import _close_exchanges

   closed = _close_exchanges(oriented)
   S = np.asarray(create_stoichiometric_matrix(closed, array_type="dense"))
   position = {r.id: i for i, r in enumerate(closed.reactions)}

   for cycle, direction in zip(cycles, directions):
       v = np.zeros(len(closed.reactions))
       for rxn_id, value in zip(cycle, direction):
           v[position[rxn_id]] = value
       assert np.abs(S @ v).max() < 1e-6

On the five *E. coli* models tested this holds exactly — the worst imbalance
across every cycle of every model was 0.

Overlapping cycles
------------------

Add ``R5``: ``B → A``. Now there are three cycles, not two, because ``R2``
and ``R3`` can run forward against ``R5`` reversed. Both enumerators find all
three in both directions; see :doc:`../functions/enumeration`.

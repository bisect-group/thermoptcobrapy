Constructing a context-specific model
=====================================

A network small enough to check by hand
---------------------------------------

Four internal reactions over four metabolites, plus two exchange reactions.
Every one of them is reversible, so the arrow is only the reference
direction — the one a positive flux runs in. ``R1``, ``R2`` and ``R3`` can
form the cycle ``A → B → C → A``, which is a TIC; ``R4`` drains ``A`` to
``D``, which the network can export:

.. code-block:: python

   from cobra import Model, Metabolite, Reaction

   model = Model("toy")          # "toy" is just the model id
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

Keeping the core without inheriting the TICs
--------------------------------------------

Ask for a CSM that keeps ``R4`` as its core:

.. code-block:: python

   >>> from thermoptcobrapy import thermoptics
   >>> extracted, blocked_core = thermoptics(model, ["R4"], tol=1e-4, seed=0)
   >>> sorted(r.id for r in extracted.reactions)
   ['EX_A', 'EX_D', 'R4']
   >>> blocked_core
   []

Three reactions survive, and ``B`` and ``C`` drop out with the TIC that was
the only thing touching them:

.. code-block:: python

   >>> sorted(m.id for m in extracted.metabolites)
   ['A', 'D']

The TIC is gone. It could have been retained — ``R1``, ``R2`` and ``R3``
are flux consistent, and an extraction requiring only flux consistency might
well have kept them — but none of them support ``R4``, and none of them can
carry thermodynamically feasible flux.

Asking for something impossible
-------------------------------

A core reaction that is thermodynamically blocked cannot be kept. Rather than
failing, it is dropped and reported:

.. code-block:: python

   >>> extracted, blocked_core = thermoptics(model, ["R1", "R4"],
   ...                                       tol=1e-4, seed=0)
   UserWarning: 1 core reaction(s) are thermodynamically blocked and have
   been dropped from the core; they are returned alongside the model.
   >>> blocked_core
   ['R1']

With ``do_tocc=False`` ThermOptCC does not run first, so nothing drops the
blocked reaction from the core and the program is infeasible:

.. code-block:: python

   >>> thermoptics(model, ["R1"], tol=1e-4, do_tocc=False, seed=0)
   ThermOptInfeasibleCoreError: None of the 1 remaining core reaction(s) can
   carry thermodynamically feasible flux ...

Only turn it off on a GEM already known to be thermodynamically flux
consistent.

What survives
-------------

Everything the source carried:

.. code-block:: python

   >>> model.reactions.R4.annotation = {"kegg.reaction": "R00001"}
   >>> extracted, _ = thermoptics(model, ["R4"], tol=1e-4, seed=0)
   >>> extracted.reactions.R4.annotation
   {'kegg.reaction': 'R00001'}

Annotations on reactions, metabolites and genes, groups, the objective and
its direction, and the order the elements were in. Gene annotations need
particular care — a gene's metadata does not survive the round trip through a
re-parsed GPR — so they are restored from the source deliberately.

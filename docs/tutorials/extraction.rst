Extracting a context-specific model
===================================

Keeping the core without inheriting the cycles
-----------------------------------------------

On the ring from :doc:`cycles`, ask for a model that keeps ``R4``:

.. code-block:: python

   >>> from thermoptcobrapy import thermoptics
   >>> extracted, blocked_core = thermoptics(model, ["R4"], tol=1e-4, seed=0)
   >>> sorted(r.id for r in extracted.reactions)
   ['EX_A', 'EX_D', 'R4']
   >>> blocked_core
   []

The ring is gone. It could have been kept — ``R1``, ``R2`` and ``R3`` are
flux consistent, and an extraction that asked only for consistency might well
have retained them — but none of it supports ``R4``, and none of it can carry
thermodynamically feasible flux.

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

With ``do_tocc=False`` the consistency check does not run first, so there is
nothing to drop the blocked reaction and the program is simply infeasible:

.. code-block:: python

   >>> thermoptics(model, ["R1"], tol=1e-4, do_tocc=False, seed=0)
   ThermOptInfeasibleCoreError: None of the 1 remaining core reaction(s) can
   carry thermodynamically feasible flux ...

Only turn it off on a model already known to be thermodynamically flux
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

On a real model
---------------

Five *E. coli* models from the protocol's Step 3, each with 40 randomly
chosen core reactions:

.. list-table::
   :header-rows: 1
   :widths: 16 18 22 22 22

   * - Model
     - Reactions
     - Extracted
     - Core kept
     - Blocked afterwards
   * - m398
     - 1523
     - 339
     - 40 / 40
     - 0
   * - m238
     - 1443
     - 461
     - 40 / 40
     - 0
   * - m463
     - 1530
     - 269
     - 40 / 40
     - 0
   * - m143
     - 1546
     - 324
     - 40 / 40
     - 0
   * - m166
     - 1202
     - 424
     - 40 / 40
     - 0

The last column is the point: each extracted model was re-checked with
:doc:`../functions/thermopt_cc` and contains no thermodynamically blocked
reactions at all.

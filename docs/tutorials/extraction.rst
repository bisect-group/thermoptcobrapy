Constructing a context-specific model
=====================================

Keeping the core without inheriting the TICs
--------------------------------------------

On the network from :doc:`cycles`, ask for a CSM that keeps ``R4``:

.. code-block:: python

   >>> from thermoptcobrapy import thermoptics
   >>> extracted, blocked_core = thermoptics(model, ["R4"], tol=1e-4, seed=0)
   >>> sorted(r.id for r in extracted.reactions)
   ['EX_A', 'EX_D', 'R4']
   >>> blocked_core
   []

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

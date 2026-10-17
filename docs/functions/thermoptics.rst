thermoptics
===========

Constructs a thermodynamically consistent context-specific model (CSM).

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermoptics

   extracted, blocked_core = thermoptics(
       model, core=["PGI", "PFK", "FBA"], tol=1e-4
   )

``blocked_core`` lists the core reactions that had to be dropped because they
are thermodynamically blocked. A warning is raised when that happens rather
than the request failing, so an imperfect core list still produces a model.

Why not extract first and check afterwards
------------------------------------------

CSM extraction normally asks for the most parsimonious sub-network carrying
non-zero flux through a given set of core reactions, and that flux is
usually required only to be *flux consistent*. This leaves a gap: a core
reaction can carry flux with support from a TIC, and the CSM then inherits
the TIC along with it.

``thermoptics`` closes the gap by requiring the flux to be thermodynamically
feasible as well, using the same indicator program as :doc:`thermopt_cc`.

How it works
------------

In rounds, through the two MILP formulations the paper calls *unBlockRxns*
and *minNetwork*. Each round

1. establishes which of the still-unplaced core reactions can carry
   thermodynamically feasible flux together, maximising their count
   (*unBlockRxns*);
2. identifies a minimal network supporting them, penalising the inclusion of
   reactions from the non-core penalty set, with the core reactions from
   step 1 required to carry non-zero flux (*minNetwork*).

Penalising non-core reactions is what yields the most parsimonious model
explaining the given biological context; without it the result could include
a large fraction of the GEM. Everything settled drops out of both the core
and the penalty set, and the next round takes what is left. A reaction joins
the network when its indicator is set *or* it carries flux, which catches
the support reactions the objective had no reason to switch on explicitly.

With ``do_tocc=True`` (the default) ThermOptCC runs first and the
thermodynamically blocked reactions leave the core. Turning it off is only
safe on a GEM already known to be thermodynamically flux consistent: a
blocked core reaction then makes the program infeasible, and
:class:`~thermoptcobrapy.ThermOptInfeasibleCoreError` is raised rather than
something wrong being returned.

``seed`` and ``time_limit``
----------------------------

.. code-block:: python

   extracted, blocked_core = thermoptics(
       model, core=["PGI", "PFK", "FBA"], tol=1e-4, seed=0, time_limit=120
   )

``seed`` makes a run reproducible, the same as in :doc:`thermopt_cc`: it
seeds both the consistency pre-check (when ``do_tocc=True``) and Fast-SNP.

``time_limit`` (default 120) is the budget given to *each* mixed-integer
solve — the two per round described above — and it is fixed, not grown the
way :doc:`thermopt_cc`'s is. A solve that times out before finding any
incumbent surfaces as
:class:`~thermoptcobrapy.ThermOptInfeasibleCoreError`, with the solver status
quoted in the message; when that status is a timeout rather than a proven
infeasibility, raising ``time_limit`` is the fix, not ``do_tocc``.

What comes back
---------------

The extracted model keeps everything the source carried — annotations on
reactions, metabolites and genes, groups, the objective, and the order the
elements were in. It is built from the kept reactions rather than copied and
cut down, which is both faster and, because a gene's metadata does not
survive the round trip through a re-parsed GPR, something that has to be done
deliberately.

On the five models tested, every CSM re-checked with zero thermodynamically
blocked reactions, which is the property the algorithm exists to provide.

Reference
---------

.. autofunction:: thermoptcobrapy.thermoptics

See also
--------

* :doc:`../tutorials/extraction` — a worked example
* :doc:`thermopt_cc` — the consistency check it runs first

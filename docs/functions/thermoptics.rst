thermoptics
===========

Builds a context-specific model that is thermodynamically feasible by
construction.

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

Why not just extract and check afterwards
-----------------------------------------

Context-specific extraction normally asks for the smallest sub-network that
keeps a given core active, and "active" is usually read as *flux consistent*.
That leaves a gap: a core reaction can be kept alive by a cycle, and the model
handed back inherits the cycle along with it.

``thermoptics`` closes the gap by asking for activity that is also
thermodynamically feasible, using the same indicator program as
:doc:`thermopt_cc`.

How it works
------------

In rounds. Each round

1. establishes which of the still-unplaced core reactions can be active
   together, maximising their count;
2. finds the cheapest set of non-core reactions supporting them, minimising
   the number switched on, with the core from step 1 forced active.

Everything settled drops out of both the core and the penalty set, and the
next round takes what is left. A reaction joins the network when its
indicator is set *or* it simply carries flux, which catches the support
reactions the objective had no reason to switch on explicitly.

With ``do_tocc=True`` (the default) the consistency check runs first and the
blocked reactions leave the core. Turning it off is only safe on a model
already known to be thermodynamically flux consistent: a blocked core
reaction then makes the program infeasible, and
:class:`~thermoptcobrapy.ThermOptInfeasibleCoreError` is raised rather than
something wrong being returned.

What comes back
---------------

The extracted model keeps everything the source carried — annotations on
reactions, metabolites and genes, groups, the objective, and the order the
elements were in. It is built from the kept reactions rather than copied and
cut down, which is both faster and, because a gene's metadata does not
survive the round trip through a re-parsed GPR, something that has to be done
deliberately.

On the five models tested, every extracted model re-checked with zero
thermodynamically blocked reactions, which is the property the whole routine
exists to provide.

Reference
---------

.. autofunction:: thermoptcobrapy.thermoptics

See also
--------

* :doc:`../tutorials/extraction` — a worked example
* :doc:`thermopt_cc` — the consistency check it runs first

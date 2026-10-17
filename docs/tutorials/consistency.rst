Detecting the blocked reactions
===============================

A network small enough to check by hand
---------------------------------------

Four internal reactions over four metabolites, plus two exchange reactions.
Every one of them is reversible, so the arrow is only the reference
direction — the one a positive flux runs in. ``R1``, ``R2`` and ``R3`` can
form the cycle ``A → B → C → A``; ``R4`` drains ``A`` to ``D``:

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

``B`` and ``C`` take part in no other reaction, so at steady state ``R1``,
``R2`` and ``R3`` must all carry the same flux — and that flux returns ``A``
to itself, producing nothing. It is a TIC. FVA predicts the maximum through
it regardless:

.. code-block:: python

   >>> from cobra.flux_analysis import flux_variability_analysis
   >>> flux_variability_analysis(model, fraction_of_optimum=0.0,
   ...                           processes=1).round(3)
         minimum  maximum
   R1    -1000.0   1000.0
   R2    -1000.0   1000.0
   R3    -1000.0   1000.0
   R4    -1000.0   1000.0
   EX_A  -1000.0   1000.0
   EX_D  -1000.0   1000.0

Running ThermOptCC
------------------

.. code-block:: python

   >>> from thermoptcobrapy import thermopt_cc
   >>> verdicts, oriented = thermopt_cc(model, tol=1e-4, seed=0)
   >>> for r in oriented.reactions:
   ...     print(f"{r.id:<6} {verdicts[r.id]}")
   R1     Blocked
   R2     Blocked
   R3     Blocked
   R4     Unblocked
   EX_A   Unblocked
   EX_D   Unblocked

``R1``, ``R2`` and ``R3`` are flux consistent under the stoichiometric and
box constraints alone — that is what FVA just reported. Adding the
thermodynamic constraints reveals them as thermodynamically blocked.

The usual alternative: loopless FVA
-----------------------------------

ThermOptCC is not the only way to reach that answer. The established
approach is loopless FVA: compute the minimum and maximum flux of every
reaction under loopless constraints, and tag as thermodynamically blocked
any reaction whose minimum and maximum are both zero. On this network it
agrees exactly:

.. code-block:: python

   >>> flux_variability_analysis(model, fraction_of_optimum=0.0,
   ...                           loopless="fastSNP", processes=1).round(3)
         minimum  maximum
   R1        0.0      0.0
   R2        0.0      0.0
   R3        0.0      0.0
   R4    -1000.0   1000.0
   EX_A  -1000.0   1000.0
   EX_D  -1000.0   1000.0

The difference is cost, and it shows on real models rather than on six
reactions. Loopless FVA solves two optimisation problems per reaction;
ThermOptCC settles many reactions per solve, dropping each one from the
search as soon as it is shown to carry flux. The paper reports ThermOptCC
faster on 95 of the 99 BiGG models tested — on average a 49-fold reduction
against Fast-SNP-based loopless FVA and roughly 27-fold against
null-space-based localised loopless FVA.

ThermOptCC also reports the stoichiometrically blocked reactions in the same
pass, and with ``get_direction=True`` the feasible direction rather than a
bare verdict. Loopless FVA gives neither without further work.

Directions
----------

.. code-block:: python

   >>> verdicts, oriented = thermopt_cc(model, tol=1e-4,
   ...                                  get_direction=True, seed=0)
   >>> verdicts["R4"]
   'Reversible'
   >>> verdicts["R1"]
   'Blocked'

Constrain ``R4`` to one direction and the verdict follows:

.. code-block:: python

   >>> model.reactions.R4.bounds = (0.0, 1000.0)
   >>> verdicts, _ = thermopt_cc(model, tol=1e-4,
   ...                           get_direction=True, seed=0)
   >>> verdicts["R4"]
   'Forward'

Use the model that comes back
-----------------------------

Reactions irreversible in the reverse direction are converted to forward
irreversible reactions before anything runs, so the verdicts are keyed to
the returned model:

.. code-block:: python

   >>> model.reactions.R1.bounds = (-1000.0, 0.0)   # lb < 0, ub = 0
   >>> verdicts, oriented = thermopt_cc(model, tol=1e-4, seed=0)
   >>> oriented.reactions.R1.bounds
   (0.0, 1000.0)
   >>> model.reactions.R1.bounds      # the input is untouched
   (-1000.0, 0.0)

The flux maps one to one onto the original's, negated. Reading a
``"Forward"`` verdict against the input model would be reading it backwards.
This returned model is the MATLAB's ``modModel``.

A forced bound can make the question unanswerable
--------------------------------------------------

Every reaction has to carry either exactly zero flux or at least ``tol``. A
reaction pinned to a growth rate can leave no room for that, and the program
becomes infeasible — which is reported as such, rather than reported as
every reaction being blocked:

.. code-block:: python

   >>> thermopt_cc(ecoli, tol=1e-4)       # biomass forced at 25.03
   ThermOptError: The consistency program is infeasible at tol=0.0001, ...

Biomass at 25.03 requires biotin synthesis to carry 6.19e-5, below the
threshold. Either lower ``tol`` — ``5e-5`` works on that model — or release
the forced bound, which is in any case appropriate if the question is which
reactions *can* carry flux rather than what happens at a fixed growth rate.
See :doc:`numerics`.

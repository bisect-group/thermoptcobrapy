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

``tol`` is the flux a reaction has to reach to count as carrying flux, and
it is the one parameter here that cannot simply be made as small as you
like; see :doc:`numerics` before choosing a value.

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

``get_direction=True`` asks a harder question, and costs accordingly. To
decide a direction, every reversible reaction is split into a forward and a
reverse half and each half is tested on its own, so the program carries
roughly twice the binaries. The split is necessary rather than wasteful: the
search drops a reaction as soon as it is shown to carry flux, so without it
a reaction proved to run forward would never be asked whether it can also
run backward. Leave it off when ``Blocked``/``Unblocked`` is all you need.

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

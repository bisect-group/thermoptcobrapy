Choosing tol, and which solver
==============================

``tol`` is the flux a reaction must reach to count as carrying flux. It is
the one parameter here that cannot simply be made as small as you like, and
the reason is worth understanding before you pick a value.

What the indicators have to mean
--------------------------------

Each reaction carries binaries :math:`z^+` and :math:`z^-`. The intended
meaning is strict:

* both off → :math:`v = 0` exactly
* :math:`z^+` on → :math:`v \ge tol`
* :math:`z^-` on → :math:`v \le -tol`

"Thermodynamically blocked" means *no assignment of these binaries lights the
reaction up*. If the binaries can lie, every verdict is worthless.

Two ways to say it
------------------

The classical encoding uses big-M coefficients, where the M is the
reaction's own flux bound:

.. math::

   v_i - ub_i z^+_i + \varepsilon z^-_i &\le 0 \\
   v_i - \underline{f}_i z^+_i - lb_i z^-_i &\ge 0

This works only if the binaries are *exactly* 0 or 1. They are not: a solver
delivers them within its integrality tolerance, and a binary sitting
:math:`\epsilon` from zero, multiplied by M, buys :math:`M\epsilon` of flux
for free. On a model with bounds of 1000 and Gurobi's floor of 1e-9, that is
1e-6 of phantom flux on every indicator — **absolute, not relative**, so it
does not shrink with the flux involved.

Observed on a real model, with the integrality tolerance at 1e-6:

.. code-block:: text

   z+ = 9.9999990e-08    z- = 0.99999990    v = 0
   row: 0 - 1000(1e-7) + 1e-4(0.9999999) = 4.46e-14  <= 0   satisfied

:math:`z^- = 1` was supposed to force :math:`v \le -10^{-4}`. The 1e-7 of
non-integrality on :math:`z^+` bought exactly :math:`10^{-4}` of slack and
cancelled it. No constraint was broken; the constraint had stopped saying
anything.

The alternative is a *native indicator constraint* — "when this binary is
set, that inequality holds" — which Gurobi and CPLEX enforce by branching.
There is no coefficient to multiply a not-quite-integral binary by, so the
slack does not arise.

Which one you get
-----------------

Chosen automatically, and you should usually leave it alone:

* **big-M** while ``tol`` clears the slack it leaves. Usually much the
  faster.
* **native indicators** when ``tol`` is finer than big-M can resolve and the
  solver offers them.

They return the same verdicts, so there is no reason to pay for indicators
until the threshold demands them. Pass ``use_indicators=True`` or ``False``
to force either.

The floors
----------

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - big-M
     - native indicators
   * - ``tol`` floor
     - ~4e-5 on bounds of 1000
     - set by ``model.tolerance``
   * - feasibility tolerance
     - must stay ≥ 1e-6
     - works at 1e-9
   * - Speed
     - faster, often much
     - slower, sometimes faster

Both are checked up front. A threshold below what the formulation can resolve
raises with the value to use, rather than being left to surface as a wrong
answer:

.. code-block:: text

   ThermOptError: tol=1e-06 is too small for this model: its largest flux
   bound is 1000, and with an integrality tolerance of 1e-09 the solver may
   grant 1e-06 of flux to a reaction whose indicator is off. ...

The bound that matters is the model's largest, so narrowing the flux bounds
of the handful of reactions left at a default of 1000 lowers the floor for
everything.

Reactions that genuinely cannot carry ``tol`` are reported blocked. That is
the honest answer at that threshold, not a failure.

Why feasTol is not set to 1e-9
------------------------------

The MATLAB sets the solver's feasibility tolerance to 1e-9. With big-M and
Gurobi 13 that produces **provably wrong optima**. On one 1543-reaction model
Gurobi proved that no reaction could be active, while a solution with 561
active reactions was feasible to 1e-12 — verified outside the solver, with
zero row violations and perfectly integral binaries.

It is a property of the big-M form, not of the solver:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Formulation
     - feasTol 1e-6
     - feasTol 1e-9
   * - big-M
     - 1367 active
     - **0 active (wrong)**
   * - native indicators
     - 1367 active
     - 1367 active

So the floor of 1e-6 applies to the big-M path only, and your model's own
tolerance is restored after the solve.

Why GLPK is rejected
--------------------

Not for being slow. GLPK returns solutions that contradict their own
indicator constraints — ``optimal``, with an indicator set on a reaction
carrying no flux. Every solution is therefore audited against the constraints
it was supposed to satisfy, and a solver that returns one it forbids raises
:class:`~thermoptcobrapy.ThermOptSolverError`:

.. code-block:: text

   ThermOptSolverError: optlang.glpk_interface returned a solution its own
   indicator constraints forbid (R1 is marked reverse but carries 0, over
   the -0.0001 required). ...

Use Gurobi or CPLEX.

A worked choice
---------------

An *E. coli* model with biomass forced at 25.03:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``tol``
     - Outcome
   * - 1e-4
     - Program infeasible — biotin synthesis must run at 6.19e-5
   * - 5e-5
     - Works. 31 of 1543 reactions blocked
   * - 1e-5
     - Works, via native indicators. 0 blocked
   * - 1e-6
     - Works, via native indicators. 0 blocked

The window exists because biomass has trace components — biotin's
stoichiometric coefficient is 2e-6 — whose supply reactions must carry flux
far below anything the growth rate itself suggests. Releasing the forced
bound avoids the squeeze entirely.

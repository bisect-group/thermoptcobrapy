thermopt_cc
===========

Detects the stoichiometrically and thermodynamically blocked reactions of a
model (and, with ``get_direction=True``, the thermodynamically feasible flux
direction of every reaction).

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_cc

   verdicts, oriented = thermopt_cc(model, tol=1e-4)
   blocked = [r for r, v in verdicts.items() if v == "Blocked"]

With ``get_direction=True`` each reaction is categorised into one of the
four groups the MATLAB reports — ``"Forward"``, ``"Reverse"``,
``"Reversible"`` or ``"Blocked"`` — instead of
``"Blocked"``/``"Unblocked"``. That costs roughly twice the work, because the
reversible reactions are split in two and each half is tested separately: the
search drops a reaction as soon as it is shown to carry flux, so without the
split a reaction proved to run forward would never be asked whether it can
also run backward.

``seed``, ``time_limit`` and ``max_time_limit``
------------------------------------------------

.. code-block:: python

   verdicts, oriented = thermopt_cc(
       model, tol=1e-4, seed=0, time_limit=120, max_time_limit=3600
   )

``seed`` makes a run reproducible. Both the consistency pre-check (SPECTRA's
own randomisation) and Fast-SNP's random projection weights use it, so the
same seed on the same model gives the same verdicts and the same basis every
time. Leave it unset and the search still gives a correct answer, just not
necessarily the same null space basis from one run to the next.

The search settles reactions in rounds, each one a mixed-integer solve over
the reactions not yet settled. A round is given ``time_limit`` seconds
(default 120); if it times out having settled nothing new, the limit grows by
60s and the same round is solved again with the longer budget. This keeps a
slow model progressing rather than failing on the first solve that needed
slightly longer.

That growth stops at ``max_time_limit`` (default 3600). If a solve is still
timing out with reactions left unsettled once the ceiling is reached,
:class:`~thermoptcobrapy.ThermOptTimeLimitError` is raised rather than
reporting the unsettled reactions as ``Blocked``. The MATLAB has no such
ceiling — it grows the limit forever — but a result that was never proven is
not a result, so this fails loudly instead of either hanging or guessing.
Raise ``max_time_limit`` on a model that needs longer, or loosen ``tol`` so
each solve is cheaper.

Why the model comes back
------------------------

The second return value is not a convenience. Before anything else runs,
every reaction irreversible in the reverse direction (``lb < 0`` and
``ub = 0``) is converted to a forward irreversible reaction: its
stoichiometric column is negated and its bounds swapped. That is a change of
variables, not a change of network, but it means "forward" refers to the
returned model rather than to the one passed in. Keying the verdicts to
anything else would be wrong. This is the MATLAB's ``modModel``.

How it works
------------

A reaction is thermodynamically consistent when some flux distribution
exists that is *both* a steady state and free of TICs. The search
maximises the number of still-unsettled reactions that can be active at once,
drops those it finds, and repeats until no more can be added.

Each reaction carries two binaries, :math:`z^+_i` and :math:`z^-_i`, meaning
"carries forward flux" and "carries reverse flux", and each reaction that can
take part in a TIC carries a continuous :math:`G_j` standing in for its Gibbs
free energy change:

.. math::

   \text{maximize}   &\quad \sum_{i \in C} (z^+_i + z^-_i) \\
   \text{subject to} &\quad S v = 0 \\
                     &\quad z^+_i = 1 \Rightarrow v_i \ge \varepsilon \\
                     &\quad z^-_i = 1 \Rightarrow v_i \le -\varepsilon \\
                     &\quad z^+_i + z^-_i \le 1 \\
                     &\quad N_{snp}^\mathsf{T} G = 0 \\
                     &\quad z^+_k = 1 \Rightarrow G_k < 0, \;
                            z^-_k = 1 \Rightarrow G_k > 0

The last three lines are the thermodynamics: a reaction proceeds in the
direction that releases energy, so :math:`G` is negative when it runs forward
and positive when it runs in reverse. :math:`G` is confined to the orthogonal
complement of the *bounded* null space, computed by Fast-SNP
(:func:`~thermoptcobrapy.fast_snp`). A TIC would require :math:`G` to
decrease around the whole cycle and return to its starting value, which that
complement does not permit — so no assignment of :math:`z` can activate a
complete TIC.

How much stronger is this than flux consistency?
------------------------------------------------

Considerably. On five *E. coli* models from the protocol's Step 3:

.. list-table::
   :header-rows: 1
   :widths: 20 20 28 32

   * - Model
     - Reactions
     - Stoichiometrically blocked
     - Thermodynamically blocked
   * - m398
     - 1523
     - 23
     - 182
   * - m238
     - 1443
     - 21
     - 170
   * - m463
     - 1530
     - 22
     - 152
   * - m143
     - 1546
     - 22
     - 179
   * - m166
     - 1202
     - 25
     - 267

Every stoichiometrically blocked reaction is also thermodynamically blocked,
as it must be. The difference — 130 to 242 reactions per model — is reactions
that FVA and flux consistency accept, because they carry flux only with
support from a TIC.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_cc

.. autofunction:: thermoptcobrapy.thermodynamically_blocked_ids

See also
--------

* :doc:`../tutorials/consistency` — worked examples
* :doc:`../tutorials/numerics` — choosing ``tol``
* :doc:`thermoptics` — extraction built on the same program

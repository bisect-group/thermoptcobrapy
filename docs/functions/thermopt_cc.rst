thermopt_cc
===========

Reports, for every reaction, which directions it can carry flux in once
cycling is excluded.

Usage
-----

.. code-block:: python

   from thermoptcobrapy import thermopt_cc

   verdicts, oriented = thermopt_cc(model, tol=1e-4)
   blocked = [r for r, v in verdicts.items() if v == "Blocked"]

With ``get_direction=True`` the verdict is one of ``"Blocked"``,
``"Forward"``, ``"Reverse"`` or ``"Reversible"`` instead of
``"Blocked"``/``"Unblocked"``. That costs roughly twice the work, because the
reversible reactions are split in two and each half is tested separately: the
search drops a reaction as soon as it is shown to carry flux, so without the
split a reaction proved to run forward would never be asked whether it can
also run backward.

Why the model comes back
------------------------

The second return value is not a convenience. Before anything else runs, every
reaction restricted to negative flux is re-oriented: its stoichiometric column
is negated and its bounds swapped. That is a change of variables, not a change
of network, but it means "forward" refers to the returned model rather than to
the one you passed in. Keying the verdicts to anything else would be wrong.
The MATLAB returns the same thing as ``modModel``.

How it works
------------

A reaction is thermodynamically consistent when some flux distribution exists
that is *both* a steady state and free of infeasible cycling. The search
maximises the number of still-unsettled reactions that can be active at once,
drops those it finds, and repeats until no more can be added.

Each reaction carries two binaries, :math:`z^+_i` and :math:`z^-_i`, meaning
"carries forward flux" and "carries reverse flux", and each cycling reaction
carries a continuous :math:`G_j` standing in for a free energy:

.. math::

   \text{maximize}   &\quad \sum_{i \in C} (z^+_i + z^-_i) \\
   \text{subject to} &\quad S v = 0 \\
                     &\quad z^+_i = 1 \Rightarrow v_i \ge \varepsilon \\
                     &\quad z^-_i = 1 \Rightarrow v_i \le -\varepsilon \\
                     &\quad z^+_i + z^-_i \le 1 \\
                     &\quad N_{snp}^\mathsf{T} G = 0 \\
                     &\quad z^+_k = 1 \Rightarrow G_k < 0, \;
                            z^-_k = 1 \Rightarrow G_k > 0

The last three lines are the thermodynamics. :math:`G` must fall when a
cycling reaction runs forward and rise when it runs backward, and it is
confined to the orthogonal complement of the *bounded* null space, computed
by Fast-SNP (:func:`~thermoptcobrapy.fast_snp`). A cycle would need :math:`G`
to decrease all the way round and return to its starting value, which that
complement does not permit — so no assignment of :math:`z` can switch a whole
cycle on.

How much stronger is this than flux consistency?
------------------------------------------------

Considerably. On five *E. coli* models from the protocol's Step 3:

.. list-table::
   :header-rows: 1
   :widths: 20 20 28 32

   * - Model
     - Reactions
     - Flux-blocked
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

Every flux-blocked reaction is also thermodynamically blocked, as it must be.
The difference — 130 to 242 reactions per model — is reactions that FVA and
plain flux consistency accept because they can carry flux, but only by
cycling.

Reference
---------

.. autofunction:: thermoptcobrapy.thermopt_cc

.. autofunction:: thermoptcobrapy.thermodynamically_blocked_ids

See also
--------

* :doc:`../tutorials/consistency` — worked examples
* :doc:`../tutorials/numerics` — choosing ``tol``
* :doc:`thermoptics` — extraction built on the same program

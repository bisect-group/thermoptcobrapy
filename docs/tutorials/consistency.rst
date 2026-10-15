Checking consistency
====================

On the ring from :doc:`cycles`
------------------------------

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

FVA gave ``R1`` a range of ±1000. The mass balance permits the ring; nothing
else does.

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

Reactions restricted to negative flux are re-oriented before anything runs,
so the verdicts are keyed to the returned model:

.. code-block:: python

   >>> model.reactions.R1.bounds = (-1000.0, 0.0)   # reverse-only
   >>> verdicts, oriented = thermopt_cc(model, tol=1e-4, seed=0)
   >>> oriented.reactions.R1.bounds
   (0.0, 1000.0)
   >>> model.reactions.R1.bounds      # the input is untouched
   (-1000.0, 0.0)

The flux maps one to one onto the original's, negated. Reading a ``"Forward"``
verdict against the input model would be reading it backwards.

A forced bound can make the question unanswerable
--------------------------------------------------

Every reaction has to be either exactly zero or at least ``tol``. A reaction
pinned to a growth rate can leave no room for that, and the program becomes
infeasible — which is reported rather than dressed up as "everything is
blocked":

.. code-block:: python

   >>> thermopt_cc(ecoli, tol=1e-4)       # biomass forced at 25.03
   ThermOptError: The consistency program is infeasible at tol=0.0001, ...

Biomass at 25.03 needs biotin synthesis running at 6.19e-5, below the
threshold. Either lower ``tol`` — ``5e-5`` works on that model — or release
the forced bound, which is anyway the right move if the question is which
reactions *can* carry flux rather than what happens at a fixed growth rate.
See :doc:`numerics`.

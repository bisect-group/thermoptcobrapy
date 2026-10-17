Installation
============

Every internal consistency check this package runs goes through
`spectra-cobra <https://github.com/bisect-group/spectra-cobra>`_, so it is a
required dependency rather than an optional one. It is not published on
PyPI, so clone and install it first:

.. code-block:: bash

   git clone https://github.com/bisect-group/spectra-cobra.git
   cd spectra-cobra
   pip install -e .
   cd ..

Then thermoptcobrapy itself:

.. code-block:: bash

   git clone https://github.com/bisect-group/thermoptcobrapy.git
   cd thermoptcobrapy
   pip install -e .

Needs Python 3.9 or newer and cobrapy 0.29 or newer.

Solvers
-------

Everything here is a mixed-integer program, and **the choice of solver is not
a question of speed**. Install Gurobi or CPLEX:

.. code-block:: bash

   pip install -e ".[gurobi]"

GLPK is rejected rather than tolerated. It has been seen returning solutions
that contradict their own indicator constraints — reporting ``optimal`` with
a reaction marked active while carrying no flux — so every solution is
checked against the constraints it was supposed to satisfy, and a solver that
returns one it forbids raises :class:`~thermoptcobrapy.ThermOptSolverError`
rather than being believed.

See :doc:`tutorials/numerics` for what that check is guarding against.

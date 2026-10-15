ThermOpt for cobrapy
====================

Flux balance analysis constrains a network with the mass balance alone, and
the mass balance has no opinion about flux that goes round in a circle. A set
of reactions whose net effect is nothing — no metabolite produced, none
consumed, nothing exchanged with the medium — can carry arbitrary flux at
steady state. Thermodynamics forbids it: the free energy would have to fall
all the way round the cycle and arrive back where it started.

These *thermodynamically infeasible cycles* are not a curiosity. They inflate
maximum fluxes, make blocked reactions look usable, and are inherited by every
context-specific model built on a network that still contains them.

``thermoptcobrapy`` is a Python port of the `ThermOptCOBRA
<https://github.com/bisect-group/Protocol_ThermOptCOBRA>`_ MATLAB package,
built on `cobrapy <https://github.com/opencobra/cobrapy>`_.

At a glance
-----------

.. code-block:: python

   from cobra.io import load_model
   from thermoptcobrapy import thermopt_cc, thermopt_enum_lp

   model = load_model("textbook")

   # Which reactions can carry thermodynamically feasible flux?
   verdicts, oriented = thermopt_cc(model, tol=1e-4)

   # Which cycles does the model contain?
   cycles, directions, cycle_rxns, oriented = thermopt_enum_lp(model)

Note that both return the model the answers refer to. Reactions restricted to
negative flux are re-oriented first, so "forward" means something different
afterwards; see :doc:`functions/thermopt_cc`.

What is in the package
----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Function
     - Does
   * - :doc:`functions/thermopt_cc`
     - Reports which directions each reaction can carry flux in, once
       cycling is excluded
   * - :doc:`functions/enumeration`
     - Lists the cycles themselves, with the relative flux through each
       reaction
   * - :doc:`functions/thermopt_flux`
     - Takes the cycling back out of a flux distribution that has it
   * - :doc:`functions/thermoptics`
     - Builds a context-specific model that is thermodynamically feasible
       by construction

Which function do I want?
-------------------------

* To find out which reactions of a model are *genuinely* usable, rather than
  usable only by cycling, use :doc:`functions/thermopt_cc`. It is strictly
  stronger than flux consistency: on the *E. coli* models tested it blocked
  130 to 242 reactions beyond what FVA-style consistency rejects.
* To see the cycles, use :doc:`functions/enumeration`. The LP enumerator is
  the cheaper; the MILP one returns tighter cycles.
* To clean up an FBA or sampling result you already have, use
  :doc:`functions/thermopt_flux`.
* To build a context-specific model without inheriting the cycles of the
  reference network, use :doc:`functions/thermoptics`.
* Before choosing ``tol``, read :doc:`tutorials/numerics`. It is the one
  parameter here that cannot simply be made small.

Citation
--------

   Thermodynamically constrained reconstruction and analysis.
   https://www.sciencedirect.com/science/article/pii/S2589004225012660

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation

.. toctree::
   :maxdepth: 2
   :caption: Functions

   functions/thermopt_cc
   functions/enumeration
   functions/thermopt_flux
   functions/thermoptics

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/cycles
   tutorials/consistency
   tutorials/extraction
   tutorials/numerics

.. toctree::
   :maxdepth: 1
   :caption: About

   api

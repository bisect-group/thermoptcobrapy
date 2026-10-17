ThermOpt for cobrapy
====================

Flux balance analysis constrains a network with the mass balance alone, and
the mass balance does not account for thermodynamically infeasible cycles
(TICs). Analogous to perpetual motion machines, these cycles violate the
second law of thermodynamics by cycling metabolites indefinitely without any
real change.

TICs lead genome-scale metabolic models (GEMs) to predict phenotypes that are
thermodynamically infeasible, and such phenotypes carry no biological
interpretation. Flux analysis methods — FBA, FVA and flux sampling alike —
often predict maximum flux through reactions involved in TICs, which
undermines a model's predictive capability through distorted flux
distributions, erroneous growth and energy predictions, unreliable gene
essentiality predictions, and compromised multi-omics integration. TICs are
also inherited by every context-specific model built on a network that still
contains them.

``thermoptcobrapy`` is a Python port of the `ThermOptCOBRA
<https://github.com/bisect-group/Protocol_ThermOptCOBRA>`_ MATLAB package,
built on `cobrapy <https://github.com/opencobra/cobrapy>`_.

At a glance
-----------

.. code-block:: python

   from cobra.io import load_model
   from thermoptcobrapy import thermopt_cc, thermopt_enum_milp

   model = load_model("textbook")

   # Thermodynamically feasible flux directions, and the blocked reactions
   verdicts, oriented = thermopt_cc(model, tol=1e-4)

   # Enumerate the TICs of the model
   cycles, directions, cycle_rxns, oriented = thermopt_enum_milp(model)

Note that both return the model the answers refer to. Reactions irreversible
in the reverse direction (``lb < 0`` and ``ub = 0``) are converted to forward
irreversible reactions first, so "forward" means something different
afterwards; see :doc:`functions/thermopt_cc`.

What is in the package
----------------------

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Function
     - Does
   * - :doc:`functions/thermopt_cc`
     - Determines the thermodynamically feasible flux directions of every
       reaction, and thereby the stoichiometrically and thermodynamically
       blocked reactions
   * - :doc:`functions/enumeration`
     - Enumerates all TICs of a model and the flux direction of the
       reactions involved
   * - :doc:`functions/thermopt_flux`
     - Detects and removes loops from a flux distribution, projecting it to
       the nearest distribution in the thermodynamically feasible flux space
   * - :doc:`functions/thermoptics`
     - Constructs a thermodynamically consistent context-specific model
       (CSM)

Which function do I want?
-------------------------

* To identify the stoichiometrically and thermodynamically blocked reactions
  of a model, use :doc:`functions/thermopt_cc`. It is strictly stronger than
  flux consistency: on the *E. coli* models tested it blocked 130 to 242
  reactions beyond what flux consistency alone rejects.
* To enumerate the TICs themselves, use :doc:`functions/enumeration`.
* To remove loops from an FBA or flux sampling result you already have, use
  :doc:`functions/thermopt_flux`.
* To construct a CSM that does not inherit the TICs of the reference GEM,
  use :doc:`functions/thermoptics`.
* Before choosing ``tol``, read :doc:`tutorials/numerics`. It is the one
  parameter here that cannot simply be made small.

Citation
--------

   Pavan Kumar S. and Nirav Pravinbhai Bhatt (2025). ThermOptCobra:
   Thermodynamically optimal construction and analysis of metabolic networks
   for reliable phenotype predictions. *iScience* **28**\ (8), 113005.
   https://doi.org/10.1016/j.isci.2025.113005

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

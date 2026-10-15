API reference
=============

A flat index of everything the package exports. The pages under *Functions*
cover the same material with explanation and worked examples; this is for
looking up a signature.

Consistency
-----------

.. autofunction:: thermoptcobrapy.thermopt_cc
   :no-index:
.. autofunction:: thermoptcobrapy.thermodynamically_blocked_ids
   :no-index:

Cycle enumeration
-----------------

.. autofunction:: thermoptcobrapy.thermopt_enum_lp
   :no-index:
.. autofunction:: thermoptcobrapy.thermopt_enum_milp
   :no-index:
.. autofunction:: thermoptcobrapy.reaction_order
   :no-index:

Flux correction
---------------

.. autofunction:: thermoptcobrapy.thermopt_flux
   :no-index:
.. autofunction:: thermoptcobrapy.cycling_removed
   :no-index:

Context-specific extraction
---------------------------

.. autofunction:: thermoptcobrapy.thermoptics
   :no-index:

Shared machinery
----------------

Useful on their own, and what the four routines above are built from.

.. autofunction:: thermoptcobrapy.fast_snp
   :no-index:
.. autofunction:: thermoptcobrapy.cycle_basis
   :no-index:
.. autofunction:: thermoptcobrapy.cycle_reaction_ids
   :no-index:
.. autofunction:: thermoptcobrapy.cycle_submodel
   :no-index:
.. autofunction:: thermoptcobrapy.prepare
   :no-index:
.. autofunction:: thermoptcobrapy.flip_reverse_only
   :no-index:
.. autofunction:: thermoptcobrapy.normalise_bounds
   :no-index:
.. autofunction:: thermoptcobrapy.convert_to_irreversible
   :no-index:
.. autofunction:: thermoptcobrapy.subset_model
   :no-index:

Errors
------

.. autoclass:: thermoptcobrapy.ThermOptError
   :no-index:
.. autoclass:: thermoptcobrapy.ThermOptSolverError
   :no-index:
.. autoclass:: thermoptcobrapy.ThermOptInfeasibleCoreError
   :no-index:
.. autoclass:: thermoptcobrapy.ThermOptTimeLimitError
   :no-index:

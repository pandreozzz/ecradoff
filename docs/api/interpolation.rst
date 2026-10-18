Interpolation tools
===================

These tools handle on-the-fly interpolation of gridded data. They are structured on serveral submodules providing different level of abstraction.


Wrapper tools
---------------------
This is the highest level API, it exposes high-level wrapper classes for 3D interpolation tasks.

.. automodule:: interp.wrapper
   :members:
   :show-inheritance:

Grids descriptions
-------------------
This includes tools to describe and handle gridded data structures.

.. automodule:: interp.grids
   :members:
   :show-inheritance:


Vertical coordinates descriptions
---------------------------------
This module provides tools to handle vertical coordinate systems and transformations.

.. automodule:: interp.verticals
   :members:
   :show-inheritance:

Interfaces to interpolation routines
------------------------------------
This is the direct interface to all Fortran interpolation routines. Also an interface to to lookup-table (LUT) interpolation is provided, though it is used by the LUT handler module (`get_lutvals()`, see :doc:`lut`)

.. automodule:: interp.interface
   :members:
   :show-inheritance:

.. Module contents
.. ---------------

.. .. automodule:: interp
..    :members:
..    :show-inheritance:
..    :undoc-members:

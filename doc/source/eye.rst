Eye submodule
====================

This module mostly being used for importing OSLO data right now.  The plan is to eventually introduce cython bindings for the C++ ray tracer so that MTFs can be generated in real time.

.. toctree::
   :maxdepth: 2

   eye_Optics.rst
   eye_navarro1999params.rst
   eye_eyeModel.rst
   
The current output of this module is a dictionary with the following key structure:

.. graphviz:: schematicEye.dot
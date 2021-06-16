Kernels
=======

Kernels are separated-out formulae for the computation of potential energy and
forces. Since they naturally appear in the same form in lots of places (for
example Lennard-Jones potential is used in almost every non-local force field),
it's appropriate to separate it out to a single class; care must be taken,
however, so as for the functions to be inlined (otherwise massive performance
costs will be incurred).

.. toctree::
   :maxdepth: 2
   :caption: Kernels

   harmonic
   lj
   stlj
   sidechainlj
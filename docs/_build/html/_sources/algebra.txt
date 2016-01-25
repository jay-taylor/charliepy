.. include:: global.rst

General Algebra
===============
In PyCox we apply the term *general algebra* quite liberally. This typically
means algebraic constructs which are not encapsulated by the theory of Coxeter
groups (for example permutations and polynomials).

.. warning::

    The implementation of the general algebra concepts are not robust. Very few
    checks are performed so that they are efficient for the task at hand. This
    means they are easily broken!

The following concepts are defined and used throughout PyCox.

.. toctree::

    permutations

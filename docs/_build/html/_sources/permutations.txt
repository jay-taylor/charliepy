.. include:: global.rst
.. module:: charliepy.algebra.permutat

Permutations
============
In PyCox a permutation is always considered to be a bijective function from the
set of integers `\{0,1,\dots,n-1\}` to itself, for some integer `n \geqslant 0`.
Such a bijective function `f` is implemented as a list of integers ``L`` such
that ``L[i]`` is simply `f(i)` for any `i \in \{0,\dots,n-1\}`. However it is
typically represented by its decomposition as a product of disjoint cycles.

.. note::

    In mathematics it is typical to consider permutations as bijective functions
    on the set `\{1,\dots,n\}`, however this is not practical in Python due to
    its indexing conventions. To obtain such functions one may simply set the
    permutation to fix 0.

Our implementation of permutations is very similar to the implementation of
permutations in `GAP`_. Indeed, many of the functions we consider below also have
analogues in GAP. However, we note the following major difference between our
permutations and those of GAP.

.. warning::

    If `f` and `g` are permutations then the composition `f\circ g` is read from
    **right to left**. In other words `(f\circ g)(i) = f(g(i))` for all `i`.
    See also ``pcx.permutat.perm.__mul__`` for more details.

Permutations are implemented through the following Python class.

.. autoclass:: Perm

Cycles
------
When dealing with permutations, two-row notation is often not very user
friendly. It can also be very inefficient when dealing with permutations
consisting of many fixed points. Thus one normally writes a permutation as a
product of disjoint cycles. The following functions allow one to produce the
cycle decomposition of a permutation and to also produce permutations from a
cycle decomposition (amongst other things).

.. autofunction:: cycledecomp
.. autofunction:: cyclestoperm
.. autofunction:: cycle
.. autofunction:: cycles
.. autofunction:: cycletype

Induced Permutations
--------------------
.. autofunction:: permonsets
.. autofunction:: restrictedperm

Miscellaneous
-------------
.. autofunction:: orderperm

.. include:: global.rst
.. currentmodule:: pycox.coxcos

:math:`\phi`-Conjugacy Classes
==============================

Let us recall that the automorphism `\phi` induces a natural equivalence
relation on the group `W`, which we call `\phi`-conjugacy. In
particular, we say two elements `x,y \in W` are `\phi`-conjugate
if there exists an element `z \in W` such that

.. math:: x = z^{-1}y\phi(z).

The equivalence classes under this relation are known as the
`\phi`-conjugacy classes of `W`, which we denote by `H^1(\phi,W)`.

.. admonition:: Definition

    We say `C \subset W\phi` is a `\phi`-conjugacy class if any
    one of the following equivalent conditions is satisfied:

    * `C` is a `W`-orbit under the natural conjugation action of `W` on
      `W\phi`.
    * `C` is a conjugacy class of `\widetilde{W}` contained in the
      coset `W\phi`.
    * the set `\{w \in W \mid (w, \phi) \in C\}` is a `\phi`-conjugacy
      class of `W`.

We now consider how we may reduce the problem of determining `\phi`-conjugacy
classes to the case where `W` is irreducible. First let us write `W` as a
direct product `W^{(1)} \times \cdots \times W^{(r)}` such that each `W^{(i)}`
is an orbit of `\phi` acting on the irreducible factors of `W`. We will denote
by `\phi^{(i)}` the restriction of `\phi` to the orbit `W^{(i)}` then the
following is obvious

.. admonition:: Lemma
    
    The natural product map 

    .. math::
        :nowrap:

        \begin{align*}
        W^{(1)} \times \cdots \times W^{(r)} &\to W\\
        (w_1,\dots,w_r) &\mapsto w_1\cdots w_r
        \end{align*}

    induces a bijection

    .. math::

        H^1(\phi^{(1)},W^{(1)}) \times \cdots H^1(\phi^{(r)},W^{(r)}) \to
        H^1(\phi, W)

In other words `W^{(i)} = W^{(i)}_1 \times \cdots
\times W^{(i)}_{m_i}` and `\phi` cyclically permutes the factors `W^{(i)}_j`.

The following function will produce the `\phi`-conjugacy classes of any
Coxeter coset.

.. autofunction:: conjugacyclasses


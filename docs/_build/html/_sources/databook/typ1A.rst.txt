Type `\mathrm{A}_n`
===================
Throughtout `(W,\mathbb{S})` will denote an irreducible Coxeter system of
type `\mathrm{A}_n`. We will also denote by `\Phi` a root system whose
associated reflection group is `W` and by `\Delta \subset \Phi^+ \subset
\Phi` a set of simple and positive roots respectively.

Root System
-----------
One may realise the indecomposabe root system of type `\mathrm{A}_n` in the
real vector space

.. math::
    V = \{(x_1,\dots,x_n) \mid x_1+\cdots+x_n = 0\} \subset \mathbb{R}^n.

In the following table we assume that `\{e_1,\dots,e_n\}` denotes the standard
basis of `\mathbb{R}^n` (i.e. `e_i = (0,\dots,0,1,\dots,0)` with 1 occurring in
the ith position).

========  ============================================================
`\Phi`    `\{e_i - e_j \mid i \neq j, 1 \leqslant i,j \leqslant n+1\}`
`\Phi^+`  `\{e_i - e_j \mid 1 \leqslant i < j \leqslant n+1\}`
`\Delta`  `\{\alpha_1 = e_1 - e_2,\dots, \alpha_n = e_n - e_{n+1}\}`
========  ============================================================

The Cartan matrix of the root system is

.. math::

    C = \begin{bmatrix}
        2 & -1 & 0 & 0 & \cdots & 0\\
        -1 & 2 & -1 & 0 & \cdots & 0\\
        0 & -1 & 2 & -1 & \cdots & 0\\
        \vdots & \vdots & \vdots & \vdots & \vdots & \vdots\\
        0 & \cdots & 0 & -1 & 2 & -1\\
        0 & \cdots & 0 & 0 & -1 & 2
        \end{bmatrix}.


Coxeter Group
-------------
In the following table we denote by `s_i` the reflection of the corresponding
simple root `\alpha_i` for `i \in \{1,\dots,n+1\}`.

=====  ====================
`W`    `\mathfrak{S}_{n+1}`
`s_i`  `(i, i+1)`
`|W|`  `(n+1)!`
=====  ====================

Here `\mathfrak{S}_{n+1}` denotes the symmetric group on `\{1,\dots,n+1\}` and
`(i,i+1) \in \mathfrak{S}_{n+1}` is a basic transposition. The Coxeter matrix of
the Coxeter system is

.. math::

    M = \begin{bmatrix}
        1 & 3 & 2 & 2 & \cdots & 2\\
        3 & 1 & 3 & 2 & \cdots & 2\\
        2 & 3 & 1 & 3 & \cdots & 2\\
        \vdots & \vdots & \vdots & \vdots & \vdots & \vdots\\
        2 & \cdots & 2 & 3 & 1 & 3\\
        2 & \cdots & 2 & 2 & 3 & 1
        \end{bmatrix}.

.. _typ-1A-coxeter-conj-classes:

Conjugacy Classes
-----------------
The cycle type of any element of the symmetric group `\mathfrak{S}_{n+1}` gives
a :term:`partition` of `n+1` which uniquely determines the conjugacy class
containing the element (see 1.2 of [JK81]_). To each partition `\lambda \vdash
n+1` we denote by `w_{\lambda}` a *very good* representative of the
corresponding conjugacy class. This may be constructed in the following way
(see 2.1 of [GM97]_).

If `\lambda = (\lambda_1,\dots,\lambda_r)` then there is a corresponding
standard parabolic subgroup `W_{\lambda}` of `W` isomorphic to
`\mathfrak{S}_{\lambda_1} \times \cdots \times \mathfrak{S}_{\lambda_r}`. The
element `w_{\lambda}` may then simply be taken as a Coxeter element of
`W_{\lambda}`. In particular we have `w_{\lambda} = w_1\cdots w_r` where

.. math::

    w_i = s_{\lambda_1+\cdots+\lambda_{i-1}+1}\cdots
    s_{\lambda_1+\cdots+\lambda_i}

Let us now write the partition `\lambda` as `(1^{a_1},\dots,n^{a_n})` for some
non-negative integers `a_i`, then we have the centraliser order is given by

.. math::
    |C_W(w_{\lambda})| = \prod_{i=1}^n i^{a_i}a_i.

(see 1.2.15 of [JK81]_).

Irreducible Characters
----------------------
There is a bijection between the set of partitions of `n+1` and the irreducible
characters of `W` which may be constructed as follows. For each
:term:`partition` `\lambda` we consider the standard parabolic subgroup
`W_{\lambda^*}` of `W` isomorphic to `\mathfrak{S}_{\lambda_1^*} \times \cdots
\times \mathfrak{S}_{\lambda_r^*}`, where `\lambda^* = (\lambda_1^*, \dots,
\lambda_r^*)` is the :term:`dual partition` of `\lambda`. The
:term:`j-induction` of the :term:`sign character` of `W_{\lambda^*}`

.. math::

    \chi_{\lambda} = j_{W_{\lambda^*}}^W(\mathrm{sgn})

is then an irreducible character of `W` (see 5.4.5 of [GP00]_). The map
`\lambda \mapsto \chi_{\lambda}` then gives the required bijection.




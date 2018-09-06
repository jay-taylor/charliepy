Type `{}^2\mathrm{A}_n`
=======================
Coxeter Coset
-------------
There is a unique non-trivial automorphism `\phi : W \to W` stabilising the set
of Coxeter generators `\mathbb{S}`. It may be realised as conjugation by the
longest element `w_0 \in W`.

`\phi`-Conjugacy Classes
------------------------
As `\phi` can be realised as an inner automorphism the `\phi`-conjugacy classes
of `W` are simply the conjugacy classes of `W`. In particular, the map `w
\mapsto ww_0` defines a bijection between the conjugacy classes and
`\phi`-conjugacy classes of `W` (see 2.9 of [GKP00]_). The map

.. math::

    \lambda \mapsto w_{\lambda} \mapsto w_{\lambda}w_0

gives a bijection between the partitions of `n+1` and the `\phi`-conjugacy
classes of `W`. Here `w_{\lambda}` is as in :ref:`typ-1A-coxeter-conj-classes`.

To construct minimal length representatives of the `\phi`-conjugacy classes we
will need to work in the symmetric group `\mathfrak{S}_{n+1}`. Let `\lambda` be
a :term:`partition` of `n+1`. We will denote by `\mu = (\mu_1,\dots,\mu_r)` a
*maximal composition* obtained from `\lambda` by rearranging the entries. In
other words, there exists `k` `(0 \leqslant k \leqslant r)` such that
`\mu_1,\dots,\mu_k` are even numbers (in any order) and `\mu_{k+1},\dots,\mu_r`
are odd numbers such that `\mu_{k+1} \geqslant \cdots \geqslant \mu_r`.

Now, to `\mu` we define an element `\sigma_{\mu} \in \mathfrak{S}_{n+1}` as
follows. Consider the list `(a_1,\dots,a_{n+1})` with `a_{2i-1} = i` and
`a_{2i} = n+1 - (i-1)` for all `1 \leqslant i \leqslant n+1` then we set

.. math::

    \sigma_{\mu} =
    \sigma_{\mu_1}\sigma_{\mu_2}\cdots\sigma_{\mu_r}

where `\sigma_{\mu_i}` is the `\mu_i`-cycle given by

.. math::

    \sigma_{\mu_i} = (a_{\mu_1+\cdots+\mu_{i-1}+1},
    a_{\mu_1+\cdots+\mu_{i-1}+2},\dots,
    a_{\mu_1+\cdots+\mu_{i-1}+\mu_i}).

By Theorem 3.3 of [GKP00]_ we have `\sigma_{\mu}w_0` is an element of
minimal length in its `\phi`-conjugacy class. From this we obtain a reduced
expression for this element by applying Algorithm A from pg.\  9 of [GP00]_.

.. warning::

    It is an open problem to determine whether these minimal length
    representatives are good in the sense of Definition 5.3 of [GKP00]_.

Irreducible Characters
----------------------
As `\phi` can be realised as an inner automorphism we have all irreducible
characters of `W` are `\phi`-stable. Thus the irreducible characters of the
coset `W\phi` are labelled by the partitions of `n+1`. Each such irreducible
character is obtained by restricting an extension `\tilde{\chi}_{\lambda}`
of `\widetilde{W} = W \rtimes \langle \phi \rangle` to the coset `W\phi`. Here
we choose Lusztig's :term:`preferred extension`. In particular, if `E` is the
simple module affording `\chi_{\lambda}` then `\phi` acts on `E` as
`(-1)^{a_{\lambda}}w_0` where `a_{\lambda}` is the a-invariant of
`\chi_{\lambda}`.




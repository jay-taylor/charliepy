*********
CharLiePy 
*********
CharLiePy (Characters of Lie type groups in Python) is a python package which
implements algorithms to carry out symbolic computations with Coxeter groups and
related mathematical structures. The ultimate goal of CharLiePy is to develop a
custom platform for carrying out computations with finite reductive groups and
their characters. It may be considered as a Python version of the original
`CHEVIE`_ package which was written in `GAP3`_ and MAPLE. This was inspired by
the Python package PyCox [Gec12]_ developed by `Meinolf Geck
<https://pnp.mathematik.uni-stuttgart.de/iaz/iaz2/geckmf/>`_. The `development
version <https://webusers.imj-prg.fr/~jean.michel/chevie/index.html>`_ of
`CHEVIE`_ is currently maintained by `Jean Michel
<https://webusers.imj-prg.fr/~jean.michel/anglais.html>`_.

CharLiePy is a mix of high level Python code and handwritten C extensions. For
instance, following `GAP`_, permutations are implemented as arrays in C to
permit faster caluclations than are possible in direct Python.  By building from
the ground up it is possible to make design choices that are specifically
tailored to working with finite reductive groups. For instance, unlike GAP,
permutations on less than 255 points are stored as 8-bit integers. This saves
memory when working with the Weyl group of :math:`\mathsf{E}_8`, which has a
natural permutation representation on 240 points.

The focus of CharLiePy is narrower than that of CHEVIE, in that we do not
anticipate working with the more general class of complex reflection groups.


Why Python?
-----------
One of the original reasons for choosing Python as the language for CharLiePy is
the success of the `Sage <http://www.sagemath.org/>`_ project. In particular,
CharLiePy may be used in Sage by simply importing CharLiePy in an interactive
session. Here are some additional reasons:

* Python is very popular with the wider scientific community (see for instance
  the `SciPy <http://www.scipy.org/>`_ project).

* It is easy to encorporate C extensions into Python, thus allowing one to
  obtain speed advantages when needed.

* As Python is a widely used modern programming language, anyone learning to use
  CharLiePy for research will simultaneously gain a useful transferable skill.

* The Python language is syntactically similar to the `GAP3`_ language, on which
  the original CHEVIE package was based. Thus it should not be difficult for
  people to switch between the two languages.

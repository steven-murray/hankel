hankel
======

.. image:: https://travis-ci.org/steven-murray/hankel.svg?branch=master
   :target: https://travis-ci.org/steven-murray/hankel
.. image:: https://coveralls.io/repos/github/steven-murray/hankel/badge.svg?branch=master
   :target: https://coveralls.io/github/steven-murray/hankel?branch=master
.. image:: https://zenodo.org/badge/19090866.svg
   :target: https://zenodo.org/badge/latestdoi/19090866
.. image:: https://img.shields.io/pypi/v/hankel.svg
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/ambv/black

Perform simple and accurate Hankel transformations using the method of
Ogata 2005.

Hankel transforms and integrals are commonplace in any area in which
Fourier Transforms are required over fields that
are radially symmetric (see
`Wikipedia <https://en.wikipedia.org/wiki/Hankel_transform>`_ for a
thorough description).
They involve integrating an arbitrary function multiplied by a Bessel
function of arbitrary order (of the first kind).
Typical integration schemes often fail because of the highly
oscillatory nature of the transform. Ogata's
quadrature method used in this package provides a fast and accurate
way of performing the integration based on
locating the zeros of the Bessel function.

Features
--------

-  Accurate and fast solutions to many Hankel integrals
-  Easy to use and re-use
-  Arbitrary order transforms
-  Built-in support for radially symmetric Fourier Transforms
-  Thoroughly tested.
-  Python 2 and 3 compatible.

Quicklinks
----------

- **Documentation:** `<https://hankel.readthedocs.io>`_
- **Quickstart+Description:** `Getting Started <https://hankel.readthedocs.io/en/latest/demos/getting_started.html>`_

Installation
------------
Either clone the repository and install locally::

    $ git clone https://github.com/steven-murray/hankel.git
    $ cd hankel/
    $ pip install -U .

Or install from PyPI:

    $ pip install hankel

The only dependencies are `numpy <https://www.numpy.org>`_,
`scipy <https://www.scipy.org>`_ and `mpmath <https://www.mpmath.org>`_.
These will be installed automatically if they are not already installed.

For instructions on testing hankel or any other development- or contribution-related
issues, see CONTRIBUTING.rst.

References
----------
Based on the algorithm provided in

    H. Ogata, A Numerical Integration Formula Based on the Bessel
    Functions, Publications of the Research Institute for Mathematical
    Sciences, vol. 41, no. 4, pp. 949-970, 2005. DOI: 10.2977/prims/1145474602

Also draws inspiration from

    Fast Edge-corrected Measurement of the Two-Point Correlation
    Function and the Power Spectrum Szapudi, Istvan; Pan, Jun; Prunet,
    Simon; Budavari, Tamas (2005) The Astrophysical Journal vol. 631 (1)
    DOI: 10.1086/496971
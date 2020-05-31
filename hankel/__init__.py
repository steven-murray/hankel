"""
General quadrature method for Hankel transformations.

Based on the algorithm provided in
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
Publications of the Research Institute for Mathematical Sciences,
vol. 41, no. 4, pp. 949-970, 2005.

.. currentmodule:: hankel

Classes
^^^^^^^

.. autosummary::
   HankelTransform
   SymmetricFourierTransform

Functions
^^^^^^^^^

.. autosummary::
   get_h

"""
from hankel.hankel import HankelTransform, SymmetricFourierTransform
from hankel.tools import get_h

try:
    from hankel._version import __version__
except ImportError:  # pragma: nocover
    # package is not installed
    pass


__all__ = ["HankelTransform", "SymmetricFourierTransform", "get_h"]

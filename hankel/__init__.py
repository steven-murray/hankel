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
from importlib.metadata import PackageNotFoundError, version

from hankel.hankel import HankelTransform, SymmetricFourierTransform
from hankel.tools import get_h

__all__ = ["HankelTransform", "SymmetricFourierTransform", "get_h"]


try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass

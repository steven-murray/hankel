Changelog
=========

dev
---
**Bugfixes**

- Fixed SymmetricFourierTransform default N to be dynamic, like HankelTransform.

v0.3.8 [08 Jan 2019]
--------------------
**Enhancements**

- Swapped to pytest from nose for all testing
- Removed warnings of overflow for cosh and sinh functions
- Improved documentation in many areas.
- Swapped out default 3.2/h for pi/h, which sounds way cooler :-)

**Bugfixes**

- Fixed an error in tests such that the wrong analytic function was being tested.
- Removed a check for deltaG in get_h, which sometimes fails spectacularly.


v0.3.7 [11 Dec 2018]
--------------------
**Bugfixes**

- Fixed issue in Py3 in which dim=2 would not run the fast function (thanks `@MuellerSeb <https://github.com/MuellerSeb>`_)

v0.3.6 [6 Aug 2018]
-------------------
**Enhancements**

- Added support for transforms with ``nu=-0.5`` (thanks `@MuellerSeb <https://github.com/MuellerSeb>`_)

v0.3.5 [8 Dec 2017]
-------------------
**Bugfixes**

- Fixed Python 3 support from v0.3.4

v0.3.4 [28 July 2017]
---------------------
**Features**

- Added get_h function to aide in determining optimal h value for a given transformation.

**Enhancements**

- Added _get_series method to quickly retrieve the summed series for the integration.
- Two updated notebook examples.

**Bugfixes**

- Moved setting of N to avoid error.

v0.3.3 [28 July 2017]
---------------------
**Features**

- Some additional tools to determine accuracy -- quick calculation of last term in sum, and evaluated range.

**Enhancements**

- Default setting of N=3.2/h, which is the maximum possible N that should be chosen, as above this, the series truncates
  due to the double-exponential convergence to the roots of the Bessel function.

**Bugfixes**

- Fixed error in cumulative sum when k is not scalar.

v0.3.2 [12 July 2017]
---------------------
**Enhancements**

- Documentation! See it at https://hankel.readthedocs.io
- Two new jupyter notebook demos (find them in the docs) by `@francispoulin <https://github.com/francispoulin>`_

**Bugfixes**

- Fixed relative import in Python 3 (tests now passing), thanks to `@louity <https://github.com/louity>`_
- Fixed docstring of SymmetricFourierTransform to have correct Fourier convention equation
- Fixed bug in choosing alternative conventions in which the fourier-dual variable was unchanged.

v0.3.1 [5 Jan 2017]
-------------------
**Bugfixes**

- Fixed normalisation for inverse transform in ``SymmetricFourierTransform``.

**Features**

- Ability to set Fourier conventions arbitrarily in ``SymmetricFourierTransform``.


v0.3.0 [4 Jan 2017]
-------------------
**Features**

- New class `SymmetricFourierTransform` which makes it incredibly easy to do arbitrary *n*-dimensional
  fourier transforms when the function is radially symmetric (includes inverse transform).
- Addition of `integrate` method to base class to perform Hankel-type integrals, which were previously
  handled by the `transform` method. This latter method is now used for actual Hankel transforms.
- Documentation!

**Enhancements**

- Addition of many tests against known integrals.
- Continuous integration
- Restructuring of package for further flexibility in the future.
- Quicker zero-finding of 1/2-order bessel functions.
- This changelog.
- Some notebooks in the devel/ directory which show how various integrals/transforms behave under
  different choices of integration steps.

---------

v0.2.2 [29 April 2016]
----------------------

**Enhancements**

- Compatibility with Python 3 (thanks to @diazona)
- Can now use with array-value functions (thanks to @diazona)

---------

v0.2.1 [18 Feb 2016]
--------------------

**Bugfixes**

- Fixed pip install by changing readme --> README

**Enhancements**

- updated docs to show dependence on mpmath

---------

v0.2.0 [10 Sep 2014]
--------------------


**Features**

* Non-integer orders supported through mpmath.

---------

v0.1.0
------
- First working version. Only integer orders (and 1/2) supported.

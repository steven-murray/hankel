Changelog
=========

Development Version
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

v0.2.2, [29 April 2016]
---------------------

**Enhancements**

- Compatibility with Python 3 (thanks to @diazona)
- Can now use with array-value functions (thanks to @diazona)

---------

v0.2.1
------
18 Feb 2016

**Bugfixes**

- Fixed pip install by changing readme --> README

**Enhancements**

- updated docs to show dependence on mpmath

---------

v0.2.0
------
10 Sep 2014

**Features**

* Non-integer orders supported through mpmath.

---------

v0.1.0
------
- First working version. Only integer orders (and 1/2) supported.

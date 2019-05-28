---
title: 'hankel: A Python library for performing simple and accurate Hankel transformations'
tags:
  - Python
  - astronomy
  - numerical analysis
authors:
  - name: Steven G. Murray
    orcid: 0000-0003-3059-3823
    affiliation: 1, 2, 3
  - name: Francis J. Poulin
    orcid: 0000-0002-7686-4089
    affiliation: 4
affiliations:
 - name: International Centre for Radio Astronomy Research (ICRAR), Curtin University,  Bentley, WA 6102, Australia
   index: 1
 - name: ARC Centre of Excellence for All-Sky Astrophysics in 3 Dimensions (ASTRO 3D)
   index: 2
 - name: School of Earth and Space Exploration, Arizona State University, Tempe, AZ, 85281, USA
   index: 3
 - name: Department of Applied Mathematics, University of Waterloo
   index: 4
date: 02 April 2019
bibliography: paper.bib
---

# Summary

The Hankel transform is a one-dimensional functional transform involving a 
Bessel-function kernel.
Moreover, it is the radial solution to an angularly symmetric Fourier 
transform of any dimension, rendering it very useful in several fields of 
application. The NASA Astronomical Data Service yields over 700 refereed 
articles including the term "hankel transform", in fields as diverse as
astronomy, geophysics, fluid mechanics, electrodynamics, thermodynamics 
and acoustics.  

As an example, in cosmology, the density field of the Universe is expected to be 
isotropic. One of the primary means of describing this field is via its Fourier 
power spectrum, or equivalently its spatial autocorrelation function. 
Due to the isotropy of the field, these can be related by an angularly symmetric 
Fourier transform, which is more simply expressed as a Hankel transform [@Szapudi2005]. 

Another example that arises in both geophysical and astrophysical contexts is in 
regards to vortices.  The radially-symmetric vortical solution to Laplace's 
equation in two, three or even higher dimensions can be performed quickly 
and accurately via the Hankel transform [@carton2001hydrodynamical].

Conceptually, computation of such problems using the Hankel transform, in 
contrast to the Fourier transform, has the advantage of reducing the problem's 
dimensionality to unity, regardless of the original dimensionality. 
Analytically, this *may* be a useful tool in solving the transform. 
Numerically, it naively promises to enhance efficiency.

Despite these advantages, the Hankel transform introduces some numerical challenges.
Most importantly, the Hankel transform is a *highly oscillatory* integral, 
especially for large values of the transformation variable, *k* (henceforth we 
will use *r* to denote the magnitude of the real-space co-ordinate). 
Highly oscillatory integrals are a topic of much interest in applied mathematics, 
and there does not exist a general optimal solution to numerically evaluate them 
double-exponential variable transformation based on the zeros of the Bessel function [@Ooura1999] 
has the property that the numerical integral converges with many fewer divisions 
compared to naively computing the transform integral. This procedure is able to 
efficiently and accurately evaluate the Hankel integral (and hence the Hankel 
in general [@huybrechs_olver_2009]. Nevertheless, [@Ogata2005] determined that a
transform) in many cases.

The purpose of ``hankel`` is to provide a dead-simple intuitive pure-Python 
interface for performing Hankel integrals and transforms, written with both 
Python 2 and 3 compatibility. It enables the accurate determination of the 
transform of a callable function at any arbitrary value of the 
transform variable, *k*, and utilises the proven method of [@Ogata2005] to do 
this efficiently. 
In addition, it recognizes that the most common application of the Hankel 
transform is in the context of the 
radial Fourier transform, and it provides an additional interface dedicated to 
making the connection between these transforms more transparent. 

The chief performance-critical components of ``hankel`` are the evaluation of the 
zeros of the Bessel function, and the sum of terms required for integration. 
The former is made efficient by a tiered approach -- using the efficient ``scipy`` 
for integer-order Bessel functions, directly returning regular arrays for 
order-1/2, and using the very accurate ``mpmath`` for all other orders. 
The latter is made efficient by utilising ``numpy``, bringing it close to C-level 
performance.  

The ``hankel`` package is thoroughly tested to ensure accuracy of transforms, 
by comparing to known analytic solutions.
These tests, supported by continuous integration, are also useful for the user 
who wishes to explore the numerical limitations of the method. Aside from 
functions which are theoretically divergent, the method can struggle to 
transform several classes of functions, including those with very sharp 
features, especially at small *r*.
The method itself has two free parameters, *h* and *N*, which respectively 
determine the resolution and upper limit of the integration grid. These can be 
modified to accurately transform any function that theoretically converges.
How to choose these values, and the estimated error of the transform under a 
given choice, are discussed in the ``hankel``'s extensive online documentation 
(and the reader is referred to [@Ogata2005] for more details).
Based on the arguments in the documentation, ``hankel`` provides an automatic, 
guided-adaptive algorithm for determination of *h* and *N*.

A particularly important limitation of ``hankel``, as currently implemented, is 
that it does *not* implement the *discrete Hankel transform*. 
That is, it provides no direct means of transforming an array of regular-spaced 
function values into "radial Fourier space" at regular-spaced *k* values. 
It is focused solely on transforming *callable* functions, so that it can 
evaluate that function at the non-regular locations required by the 
double-exponential transform of [@Ogata2005].
Extensions to discrete Hankel transforms (and even *fast* Hankel transforms) are
 envisioned for v2.0 of ``hankel``.
 

# Acknowledgements

The authors acknowledge helpful contributions from Sebastian Mueller during the 
construction of this code. Parts of this research were supported by the Australian 
Research Council Centre of Excellence for All Sky Astrophysics in 3 Dimensions 
(ASTRO 3D), through project number CE170100013. 
FJP would like to thank NSERC for research funding during the time of this research. 

# References
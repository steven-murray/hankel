------
hankel
------

Perform simple and accurate Hankel transformations using the method of Ogata 2005.

Hankel transforms and integrals are commonplace in any area in which Fourier Transforms
are required over fields that are radially symmetric (see `Wikipedia <https://en.wikipedia.org/wiki/Hankel_transform>`_
for a thorough description). They involve integrating an arbitrary function
multiplied by a Bessel function of arbitrary order (of the first kind).
Typical integration schemes often fall over
because of the highly oscillatory nature of the transform.
Ogata's quadrature method used in this package
provides a fast and accurate way of performing the integration based on locating
the zeros of the Bessel function.

Installation
------------
Either clone the repository at github.com/steven-murray/hankel and use
``python setup.py install``, or simply install using ``pip install hankel``.

The only dependencies are numpy, scipy and mpmath (as of v0.2.0).

Usage
-----

Setup
+++++
This implementation is set up to allow efficient calculation of multiple
functions *f(x)*. To do this, the format is class-based, with the main object
taking as arguments the order of the Bessel function, and the number and size
of the integration steps (see Limitations_ for discussion about how to choose
these key parameters).

For any general integration or transform of a function, we perform the following
setup:

.. code:: python

    from hankel import HankelTransform     # Import the basic class

    ht = HankelTransform(nu= 0,            # The order of the bessel function
                         N = 120,          # Number of steps in the integration
                         h = 0.03)         # Proxy for "size" of steps in integration

Alternatively, each of the parameters has defaults, so you needn't pass any.
The order of the bessel function will be defined by the problem at hand, while the other
arguments typically require some exploration to set them optimally.


Integration
+++++++++++
A Hankel-type integral is the integral

.. image:: https://rawgit.com/steven-murray/hankel/None/svgs/de0f311f07c2e7b1634b56eef090e0b2.svg?invert_in_darkmode
    :align: middle
    :width: 126.78105pt
    :height: 38.2239pt

Having set up our transform with `nu = 0`, we may wish to perform this integral for *f(x) = 1*.
To do this, we do the following:

.. code:: python

    f = lambda x : 1   # Create a function which identically 1.
    ht.integrate(f)    # Should give (1.0000000000003544, -9.8381428368537518e-15)


The correct answer is 1, so we have done quite well. The second element of the 
returned result is an estimate of the error (it is the last term in the
summation). The error estimate can be omitted using the argument
``ret_err=False``.

We may now wish to integrate a different function, say <img src="https://rawgit.com/steven-murray/hankel/None/svgs/ab2aba38e3df024ffc0b7418ccede184.svg?invert_in_darkmode" align=middle width=75.252375pt height=26.70657pt/>. We can do this
directly with the same object, without re-instantiating (avoiding unnecessary recalculation):

.. code:: python

    f = lambda x : x/(x**2 + 1)
    ht.integrate(f)               # Should give (0.42098875721567186, -2.6150757700135774e-17)

The analytic answer here is <img src="https://rawgit.com/steven-murray/hankel/None/svgs/19f9b8feb2c5698a0de2d8df238db319.svg?invert_in_darkmode" align=middle width=109.57122pt height=24.56553pt/>. The accuracy could be increased by
creating `ht` with a higher number of steps `N`, and lower stepsize `h` (see Limitations_).

Transforms
++++++++++
The Hankel transform is defined as

<p align="center"><img src="https://rawgit.com/steven-murray/hankel/None/svgs/c26683fe1415e9a3ac048789dbb84a86.svg?invert_in_darkmode" align=middle width=194.67855pt height=38.2239pt/></p>


We see that the Hankel-type integral is the Hankel transform of <img src="https://rawgit.com/steven-murray/hankel/None/svgs/a5c3ca476cddee18a9f95b93ee03a74e.svg?invert_in_darkmode" align=middle width=46.40229pt height=24.56553pt/> with <img src="https://rawgit.com/steven-murray/hankel/None/svgs/7eb22be4bf74527b54b6d60938478147.svg?invert_in_darkmode" align=middle width=39.101865pt height=22.74591pt/>.
To perform this more general transform, we must supply the <img src="https://rawgit.com/steven-murray/hankel/None/svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.041505pt height=22.74591pt/> values. Again, let's
use our previous function, <img src="https://rawgit.com/steven-murray/hankel/None/svgs/ab2aba38e3df024ffc0b7418ccede184.svg?invert_in_darkmode" align=middle width=75.252375pt height=26.70657pt/>:

.. code:: python

    import numpy as np              # Import numpy
    k = np.logspace(-1,1,50)        # Create a log-spaced array of k from 0.1 to 10.
    ht.transform(f,k,ret_err=False) # Return the transform of f at k.

Fourier Transforms
++++++++++++++++++
One of the most common applications of the Hankel transform is to solve the `radially symmetric
*n*-dimensional Fourier transform <https://en.wikipedia.org/wiki/Hankel_transform#Relation_to_the_Fourier_transform_.28radially_symmetric_case_in_n-dimensions.29>`_:

<p align="center"><img src="https://rawgit.com/steven-murray/hankel/None/svgs/a406d349d3c2e83a25072dc0dd3c3e52.svg?invert_in_darkmode" align=middle width=332.3529pt height=40.66128pt/></p>

We provide a specific class to do this transform, which takes into account the various normalisations and substitutions
required, and also provides the inverse transform. The procedure is similar to the basic `HankelTransform`, but
we provide the number of dimensions, rather than the Bessel order directly. Say we wish to find the Fourier transform
of <img src="https://rawgit.com/steven-murray/hankel/None/svgs/3967151b686bd92ca565ce7ca2ef9f14.svg?invert_in_darkmode" align=middle width=76.46067pt height=24.56553pt/> in 3 dimensions:

.. code:: python

    from hankel import SymmetricFourierTransform
    ft = SymmetricFourierTransform(ndim=3, N = 200, h = 0.03)

    f = lambda r : 1./r
    ft.transform(f,k, ret_err=False)

To do the inverse transformation (which is different by a normalisation constant), merely supply `inverse=True` to the
`.transform()` method.


Limitations
-----------
Efficiency
++++++++++
An implementation-specific limitation is that the method is not perfectly
efficient in all cases. Care has been taken to make it efficient in the general 
sense. However, for specific orders and functions, simplifications may be made
which reduce the number of trigonometric functions evaluated. For instance,
for a zeroth-order spherical transform, the weights are analytically always identically
1. 

Lower-Bound Convergence
+++++++++++++++++++++++
In terms of limitations of the method, they are very dependent on the form of the
function chosen. Notably, functions which tend to infinity at x=0 will be poorly
approximated in this method, and will be highly dependent on the step-size
parameter, as the information at low-x will be lost between 0 and the first step.
As an example consider the simple function <img src="https://rawgit.com/steven-murray/hankel/None/svgs/c67b8cebef81aeb9e88f320fe5a22d6d.svg?invert_in_darkmode" align=middle width=93.225495pt height=24.99552pt/> with a 1/2 order bessel function.
The total integrand tends to 1 at x=0, rather than 0:

.. code:: python

    f = lambda x: 1/np.sqrt(x)
    h = HankelTransform(0.5,120,0.03)
    h.integrate(f)  #(1.2336282286725169, 9.1467916948046785e-17)

The true answer is <img src="https://rawgit.com/steven-murray/hankel/None/svgs/3547c67a9cab7f8c7104501d28de4874.svg?invert_in_darkmode" align=middle width=723.97545pt height=521.06472pt/>f(x) = x^{1/2}<img src="https://rawgit.com/steven-murray/hankel/None/svgs/c2749a6ad216bc00145da6f496de7a25.svg?invert_in_darkmode" align=middle width=319.159995pt height=24.66849pt/>f(x) = x^{0.4}$ as an example of a slowly converging function, and use our "hi-res"
setup from the previous section:

.. code:: python

    h = HankelTransform(0.5,700,0.001)
    f = lambda x : x**0.4
    h.integrate(f)   # (0.53678277933471386, -1.0590954621246349)

The analytic result is 0.8421449 -- very far from our result. Note that in this case, the error estimate itself is a
good indication that we haven't reached convergence. We could try increasing `N`:

.. code:: python

    h = HankelTransform(0.5,10000,0.001)
    h.integrate(f,ret_err=False)/0.8421449 -1     ## 7.128e-07

This is very accurate, but quite slow. Alternatively, we could try increasing `h`:

.. code:: python

    h = HankelTransform(0.5,700,0.03)
    h.integrate(f,ret_err=False)/0.8421449 -1     ## 0.00045616

Not quite as accurate, but still far better than a percent for a hundredth of the cost!

There are some notebooks in the devel/ directory which toy with some known integrals, and show how accurate different
choices of `N` and `h` are. They are interesting to view to see some of the patterns.


References
----------
Based on the algorithm provided in 

   H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
   Publications of the Research Institute for Mathematical Sciences, 
   vol. 41, no. 4, pp. 949-970, 2005.

Also draws inspiration from 

   Fast Edge-corrected Measurement of the Two-Point Correlation Function and the Power Spectrum
   Szapudi,  Istvan;  Pan,  Jun;  Prunet,  Simon;  Budavari,  Tamas (2005)
   The Astrophysical Journal	vol. 631 (1)

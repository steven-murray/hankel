------
hankel
------

Perform simple and accurate Hankel transformations using the method of Ogata 2005.

Hankel transforms and integrals are commonplace in any area in which Fourier Transforms
are required over fields that are radially symmetric (see (Wikipedia)[https://en.wikipedia.org/wiki/Hankel_transform]
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

```python
from hankel import HankelTransform     # Import the basic class

ht = HankelTransform(nu= 0,            # The order of the bessel function
                     N = 120,          # Number of steps in the integration
                     h = 0.03)         # Proxy for "size" of steps in integration
```

Alternatively, each of the parameters has defaults, so you needn't pass any.
The order of the bessel function will be defined by the problem at hand, while the other
arguments typically require some exploration to set them optimally.


Integration
+++++++++++
A Hankel-type integral is the integral

$$ \int_0^\infty f(x) J_\nu(x) dx. $$

Having set up our transform with `nu = 0`, we may wish to perform this integral for *f(x) = 1*.
To do this, we do the following:

```python
f = lambda x : 1   # Create a function which identically 1.
ht.integrate(f)    # Should give (1.0000000000003544, -9.8381428368537518e-15)
```

The correct answer is 1, so we have done quite well. The second element of the 
returned result is an estimate of the error (it is the last term in the
summation). The error estimate can be omitted using the argument
``ret_err=False``.

We may now wish to integrate a different function, say $x/(x^2 + 1)$. We can do this
directly with the same object, without re-instantiating (avoiding unnecessary recalculation):

```python
f = lambda x : x/(x**2 + 1)
ht.integrate(f)               # Should give (0.42098875721567186, -2.6150757700135774e-17)
```

The analytic answer here is $K_0(1) = 0.4210$. The accuracy could be increased by
creating `ht` with a higher number of steps `N`, and lower stepsize `h` (see Limitations_).

Transforms
++++++++++
The Hankel transform is defined as

$$ F(k) = \int_0^\infty f(r) J_\nu(kr) r dr. $$

We see that the Hankel-type integral is the Hankel transform of $f(r)/r$ with $k=1$.
To perform this more general transform, we must supply the $k$ values. Again, let's
use our previous function, $x/(x^2 + 1)$:

```python
import numpy as np              # Import numpy
k = np.logspace(-1,1,50)        # Create a log-spaced array of k from 0.1 to 10.
ht.transform(f,k,ret_err=False) # Return the transform of f at k.
```

Fourier Transforms
++++++++++++++++++
One of the most common applications of the Hankel transform is to solve the (radially symmetric
*n*-dimensional Fourier transform)[https://en.wikipedia.org/wiki/Hankel_transform#Relation_to_the_Fourier_transform_.28radially_symmetric_case_in_n-dimensions.29]:

$$ F(k) = \frac{(2\pi)^{n/2}}{k^{n/2-1}} \int_0^\infty r^{n/2-1} f(r) J_{n/2-1}(kr)r dr. $$

We provide a specific class to do this transform, which takes into account the various normalisations and substitutions
required, and also provides the inverse transform. The procedure is similar to the basic `HankelTransform`, but
we provide the number of dimensions, rather than the Bessel order directly. Say we wish to find the Fourier transform
of $f(r) = 1/r$ in 3 dimensions:

```python
from hankel import SymmetricFourierTransform
ft = SymmetricFourierTransform(ndim=3, N = 200, h = 0.03)

f = lambda r : 1./r
ft.transform(f,k, ret_err=False)
```

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
As an example consider the simple function $f(x) = 1/\sqrt{x}$ with a 1/2 order bessel function.
The total integrand tends to 1 at x=0, rather than 0:

```python
f = lambda x: 1/np.sqrt(x)
h = HankelTransform(0.5,120,0.03)
h.integrate(f)  #(1.2336282286725169, 9.1467916948046785e-17)
```

The true answer is $\sqrt{pi/2}, which is a difference of about 1.6%. Modifying the step
size and number of steps to gain accuracy we find::

```python
h = HankelTransform(0.5,700,0.001)
h.integrate(f)   #(1.2523045156429067, -0.0012281146007910256)
```
This has much better than percent accuracy, but uses 5 times the amount
of steps. The key here is the reduction of h to "get inside" the low-x information.
This limitation is amplified for cases where the function really does tend to
infinity at x=0, rather than a finite positive number, such as f(x) = 1/x.
Clearly the integral becomes non-convergent for some *f(x)*, in which case
the numerical approximation can never be correct.

Upper-Bound Convergence
+++++++++++++++++++++++
If the function *f(x)* is monotonically increasing, or at least very slowly decreasing, then higher and higher zeros
of the Bessel function will be required to capture the convergence. Often, it will be the case that if this is so, the
amplitude of the function is low at low *x*, so that the step-size `h` can be increased to facilitate this. Otherwise,
the number of steps `N` can be increased.

For example, the 1/2-order integral supports functions that are increasing up to $f(x) = x^{1/2}$ and no more
(otherwise they diverge). Let's use $f(x) = x^{0.4}$ as an example of a slowly converging function, and use our "hi-res"
setup from the previous section:

```python
h = HankelTransform(0.5,700,0.001)
f = lambda x : x**0.4
h.integrate(f)   # (0.53678277933471386, -1.0590954621246349)
```

The analytic result is 0.8421449 -- very far from our result. Note that in this case, the error estimate itself is a
good indication that we haven't reached convergence. We could try increasing `N`:

```python
h = HankelTransform(0.5,10000,0.001)
h.integrate(f,ret_err=False)/0.8421449 -1     ## 7.128e-07
```

This is very accurate, but quite slow. Alternatively, we could try increasing `h`:

```python
h = HankelTransform(0.5,700,0.03)
h.integrate(f,ret_err=False)/0.8421449 -1     ## 0.00045616
```

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

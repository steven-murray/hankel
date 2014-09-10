------
hankel
------

This is a very simple module designed to implement the Hankel transformation
following the method of Ogata 2005. 

It is a fast and accurate way of integrating functions of the form f(x)J(x),
where f is an arbitrary slowly decreasing function, and J(x) is an arbitrary 
Bessel function of the first kind.

Installation
------------
Either clone the repository at github.com/steven-murray/hankel and use
``python setup.py install``, or simply install using ``pip install hankel``.

The only dependencies are numpy and scipy.

Usage
-----
This implementation is set up to allow efficient calculation of multiple
functions f(x). To do this, the format is class-based, with the main object 
taking as arguments the order of the Bessel function, and the number and size
of the integration steps. For example, to integrate the function J_0(x) (ie.
f(x) = 1, cf. Ogata's paper) one would do the following::
   
   from hankel import HankelTransform
   f = lambda x: 1  #Define the input function f(x)
   h = HankelTransform(nu=0,N=120,h=0.03)  #Create the HankelTransform instance
   h.transform(f)  #Should give [1.0000000000003544, -9.8381428368537518e-15]
   
The correct answer is 1, so we have done quite well. The second element of the 
returned result is an estimate of the error (it is the last term in the
summation). Here we used 120 steps of size 0.03. Different applications will
need to tune these parameters to get best results. In the above example, one
may modify the function ``f`` and re-call ``h.transform(f)`` without re-instantiating.
This avoids unnecessary recalculation.

Also included in the module is a subclass called ``SphericalHankelTransform``.
This is dedicated to integrating functions of the form f(x)j(x), where j(x) is 
a *spherical* Bessel function of arbitrary order. It is called in exactly the
same way. An example::

	from hankel import SphericalHankelTransform
	f = lambda x: x/(x**3+1)  #Define the input function f(x)
   	h = SphericalHankelTransform(nu=0,N=500,h=0.005)  #Create the HankelTransform instance
   	h.transform(f)  #Should give [0.61092293340214776, -1.4163951324130801e-14]
   	
Mathematica gives the answer as 0.610913.

Limitations
-----------
An implementation-specific limitation is that the method is not perfectly
efficient in all cases. Care has been taken to make it efficient in the general 
sense. However, for specific orders and functions, simplifications may be made
which reduce the number of trigonometric functions evaluated. For instance,
for a zeroth-order spherical transform, the weights are analytically always identically
1. 

In terms of limitations of the method, they are very dependent on the form of the
function chosen. Notably, functions which tend to infinity at x=0 will be poorly
approximated in this method, and will be highly dependent on the step-size
parameter, as the information at low-x will be lost between 0 and the first step.
As an example consider the simple function f(x) = 1 with a zeroth order spherical
bessel function. This tends to 1 at x=0, rather than 0:: 

   f = lambda x: 1
   h = SphericalHankelTransform(0,120,0.03)
   h.transform(f) #[1.5461236955707951, -3.5905712375161296e-16] 
   
The true answer is pi/2, which is a difference of about 3%. Modifying the step
size and number of steps to gain accuracy we find::

   h = SphericalHankelTransform(0,10000,0.0001)
   h.transform(f) #[1.5706713512229455, -0.00010492204442285768]   
   
This has much better than percent accuracy, but uses almost 100 times the amount
of steps. The key here is the reduction of h to "get inside" the low-x information.
This limitation is amplified for cases where the function really does tend to
infinity at x=0, rather than a finite positive number, such as f(x) = 1/x.

History
-------
v0.2.0 -- 10 Sep 2014. 
		  Non-integer orders supported through mpmath.
		  
v0.1.0 -- First working version. Only integer orders (and 1/2) supported.

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

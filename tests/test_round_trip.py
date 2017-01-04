'''
This module provides tests that ensure that the forward then backward transform ends up with the original.

For simplicity and because of common use (for me at least), we do this with a broken powerlaw, which turns
over at some x0. This *can* ensure convergence of the integral, as long as the power-law exponents are chosen
correctly (some here are commented out because they don't fit the convergence criteria at high x).

There is a corresponding notebook in devel/ that runs each of these functions through a grid of N and h,
showing the pattern of accuracy. This could be useful for finding the correct numbers to choose for these
for other unknown functions.
'''

import numpy as np
import inspect
import os
LOCATION = "/".join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))).split("/")[:-1])
import sys
sys.path.insert(0, LOCATION)

from hankel import SymmetricFourierTransform

from scipy.interpolate import InterpolatedUnivariateSpline as spline


class TestRoundTrip(object):

    def do_roundtrip_broken_pl(self,s,t,x0,ndim,N,h,r):
        f = lambda x : x**s/(x**t + x0)
        ht = SymmetricFourierTransform(ndim=ndim,N=N,h=h)

        k = np.logspace(np.log10(ht.x.min()/r.max()/10.0), np.log10(10*ht.x.max()/r.min()), 1000)
        resk = ht.transform(f, k, False)
        spl = spline(k, resk)
        res = ht.transform(spl, r, False, inverse=True)

        print("Analytic: ", f(r), ", Numerical: ", res)
        assert np.allclose(res,f(r),atol=1e-2)

    def test_roundtrip_broken_pl(self):
        r = np.array([0.1,1,10.])
        trials = [[0, 1, 1, 2, 1100, 10**-2.5, r],
                  [0, 2, 1, 2, 800,10**-2.6, r],
                  [0, 3, 1, 2, 1500, 10**-3.5, r[:2]],     # Doesn't work for high k!
                  # [1, 1, 1, 2, 4000, 10**-3, r],           # Doesn't seem to work for any k. SHOULDN'T CONVERGE
                  [1, 2, 1, 2, 250, 10**-2, r],
                  [1, 3, 1, 2, 1200, 0.0008, r],
                  # [2, 1, 1, 2, 5000, 10**-3, r],           # Never works. SHOULDN"T CONVERGE
                  # [2, 2, 1, 2, 5000, 10**-3, r],           # Never works. SHOULDN'T CONVERGE
                  [2, 3, 1, 2, 1200, 10**-2.5, r],
                  # [3, 1, 1, 2, 5000, 10**-3, r],           # Never works. SHOULDN'T CONVERGE
                  # [3, 2, 1, 2, 5000, 10**-3, r],           # Never works. SHOULDN'T CONVERGE
                  # [3, 3, 1, 2, 5000, 10**-3, r],           # Never works. SHOULDN'T CONVERGE
#
                  [0, 1, 10, 2, 1200, 10**-2.6, r[1:]],    # Doesn't work for low k
                  [0, 2, 10, 2, 250,  10**-2.0, r[1:]],    # Doesn't work for low k
                  [0, 3, 10, 2, 1000, 10**-3, r],
                  # [1, 1, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
                  [1, 2, 10, 2, 1100, 10**-2.5, r],
                  [1, 3, 10, 2, 250, 10**-2.2, r[1:]],     # Doesn't work for low k
                  # [2, 1, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
                  # [2, 2, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
                  [2, 3, 10, 2, 700, 10**-2.3, r[1:]],     # Doesn't work for low k
                  # [3, 1, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
                  # [3, 2, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
                  # [3, 3, 10, 2, 5000, 10**-3, r],          # Doesn't work. SHOULDN'T CONVERGE
#
                  # [0, 1, 1, 3, 5000, 10**-3, r],           # Out of range. SHOULDN'T CONVERGE
                  [0, 2, 1, 3, 700, 10**-2.4, r],
                  [0, 3, 1, 3, 250, 10**-2, r],
                  # [1, 1, 1, 3, 5000, 10**-3, r],           # Out of range. SHOULDN'T CONVERGE
                  # [1, 2, 1, 3, 5000, 10**-3, r],           # Out of range. SHOULDN'T CONVERGE
                  [1, 3, 1, 3, 2000, 10**-2.9, r],
                  # [2, 1, 1, 3, 5000, 10**-3, r],           # Out of range. SHOULDN'T CONVERGE
                  # [2, 2, 1, 3, 5000, 10**-3, r],           # Out of range. SHOULDN'T CONVERGE
                  # [2, 3, 1, 3, 50, 0.05, r],               # Out of range. SHOULDN'T CONVERGE
                  # [3, 1, 1, 3, 50, 0.05, r],               # Out of range. SHOULDN'T CONVERGE
                  # [3, 2, 1, 3, 50, 0.05, r],               # Out of range. SHOULDN'T CONVERGE
                  # [3, 3, 1, 3, 50, 0.05, r],               # Out of range. SHOULDN'T CONVERGE
                  #
                  # [0, 1, 10, 3, 50, 0.05, r],
                  [0, 2, 10, 3, 2200, 10**-2.8, r],
                  [0, 3, 10, 3, 400, 10**-2.4, r],
                  # [1, 1, 10, 3, 50, 0.05, r],
                  # [1, 2, 10, 3, 50, 0.05, r],
                  [1, 3, 10, 3, 5000, 10**-3, r],                # Bad at low k
                  # [2, 1, 10, 3, 50, 0.05, r],
                  # [2, 2, 10, 3, 50, 0.05, r],
                  # [2, 3, 10, 3, 50, 0.05, r],
                  # [3, 1, 10, 3, 50, 0.05, r],
                  # [3, 2, 10, 3, 50, 0.05, r],
                  # [3, 3, 10, 3, 50, 0.05, r],

                  ]

        for s,t,x0,ndim,N,h,r in trials:
            yield self.do_roundtrip_broken_pl, s,t,x0,ndim,N,h,r
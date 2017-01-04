'''
This module provides some tests of the transform, based on analytic transforms provided by Wikipedia.

Note that most of the tests in here are mirrored in test_known_integrals, in which k is set to 1.

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



from scipy.special import k0, gamma, gammainc, gammaincc

gammainc_ = lambda a,x : gamma(a)*gammainc(a,x)
gammaincc_ = lambda a,x : gamma(a)*gammaincc(a,x)

from hankel import HankelTransform


class TestAnalyticTransforms(object):

    def powerlaw(self,s,nu,k,N,h):
        """
        Test f(r) = 1/r, nu=0
        """

        ht = HankelTransform(nu=nu, N=N, h=h)
        ans = ht.transform(lambda x: x**s, k, False, False)
        if nu-s <= 0 and (nu-s)%2==0:
            raise Exception("Can't have a negative integer for gamma")

        anl = 2**(s+1) * gamma(0.5*(2+nu+s))/k**(s+2)/gamma(0.5*(nu-s))

        print("Numerical Result: ", ans, " (required %s)"%anl)
        assert np.isclose(ans,anl,rtol=1e-3)

    def test_powerlaw(self):
        trials = [[-2, 1, 0.01, 300, 10**-3.2],
                  [-2, 1, 1, 300, 10**-3.2],
                  [-2, 1, 10.0, 300, 10**-3.2],
                  [-2, 2, 0.01, 200, 10**-2.],
                  [-2, 2, 1, 200, 10**-2.],
                  [-2, 2, 10.0, 200, 10**-2.],
                  [-1, 0, 0.01, 50, 0.05],
                  [-1, 0, 1, 50, 0.05],
                  [-1, 0, 10.0, 50, 0.05],
                  [-1, 1, 0.01, 50, 0.05],
                  [-1, 1, 1, 50, 0.05],
                  [-1, 1, 10.0, 50, 0.05],
                  [-1, 2, 0.01, 50, 0.05],
                  [-1, 2, 1, 50, 0.05],
                  [-1, 2, 10.0, 50, 0.05],
                  [1, 0, 0.01, 150, 10**-1.5],
                  [1, 0, 1, 150, 10**-1.5],
                  [1, 0, 10.0, 150, 10**-1.5],
                  [1, 2, 0.01, 50, 0.05],
                  [1, 2, 1, 50, 0.05],
                  [1, 2, 10.0, 50, 0.05],
                  [2, 1, 0.01, 100, 10**-1.5],
                  [2, 1, 1, 100, 10**-1.5],
                  [2, 1, 10.0, 100, 10**-1.5],
                  ]
        for s,nu,k,N,h in trials:
            yield self.powerlaw, s, nu, k, N,h
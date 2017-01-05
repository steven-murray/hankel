'''
This module provides some tests for the Symmetric Fourier Transforms.

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

from hankel import SymmetricFourierTransform


def powerlaw(s,ndim,k,N,h):
    ht = SymmetricFourierTransform(ndim=ndim, N=N, h=h)
    ans = ht.transform(lambda x: x**s, k, False, False)

    nu = ndim/2. - 1
    s += nu
    if nu-s <= 0 and (nu-s)%2==0:
        raise Exception("Can't have a negative integer for gamma")

    anl = (2*np.pi)**(ndim/2.) * 2**(s+1) * gamma(0.5*(2+nu+s))/k**(s+2)/gamma(0.5*(nu-s))/ k**nu

    print("Numerical Result: ", ans, " (required %s)"%anl)
    assert np.isclose(ans,anl,rtol=1e-3)

def test_powerlaw():
    trials = [#[-2, 4, 0.01, 300, 10 ** -3.2],
    #           [-2, 4, 1, 300, 10 ** -3.2],
    #           [-2, 4, 10.0, 300, 10 ** -3.2],
    #           [-2, 6, 0.01, 200, 10 ** -2.],
    #           [-2, 6, 1, 200, 10 ** -2.],
    #           [-2, 6, 10.0, 200, 10 ** -2.],
              [-1, 2, 0.01, 50, 0.05],
              [-1, 2, 1, 50, 0.05],
              [-1, 2, 10.0, 50, 0.05],
              [-1, 4, 0.01, 50, 0.05],
              [-1, 4, 1, 50, 0.05],
              [-1, 4, 10.0, 50, 0.05],
              [-1, 6, 0.01, 50, 0.05],
              [-1, 6, 1, 50, 0.05],
              [-1, 6, 10.0, 50, 0.05],
              [1, 2, 0.01, 150, 10 ** -1.5],
              [1, 2, 1, 150, 10 ** -1.5],
              [1, 2, 10.0, 150, 10 ** -1.5],
              [1, 3, 0.01, 150, 10 ** -1.5],
              [1, 3, 1, 150, 10 ** -1.5],
              [1, 3, 10.0, 150, 10 ** -1.5],
              [-1, 3, 0.01, 150, 10 ** -1.5],
              [-1, 3, 1, 150, 10 ** -1.5],
              [-1, 3, 10.0, 150, 10 ** -1.5],

        # [1, 6, 0.01, 50, 0.05],
              # [1, 6, 1, 50, 0.05],
              # [1, 6, 10.0, 50, 0.05],
              ]
    for s, nu, k, N, h in trials:
        yield powerlaw, s, nu, k, N, h
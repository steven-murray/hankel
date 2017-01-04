'''
This module provides some tests of the integrator to known integrals.

Note, these are not the transformations, just the plain integrals,  :math:`\int_0^\infty f(x) J_\nu(x) dx`

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


class TestAnalyticIntegrals_Order0(object):

    def test_f_unity(self):
        """
        Test f(x) = 1, nu=0

        This test is done in the Ogata (2005) paper, section 5"
        """

        ht = HankelTransform(nu=0, N=50, h=0.1)
        ans = ht.integrate(lambda x: 1, False, False)
        print "Numerical Result: ", ans, " (required %s)"%1
        assert np.isclose(ans,1,rtol=1e-3)

    def test_f_x_on_x2(self):
        """
        Test f(x) = x/(x**2 + 1), nu=0

        This test is done in the Ogata (2005) paper, section 5"
        """
        ht = HankelTransform(nu=0, N=50, h=10**-1.5)


        ans = ht.integrate(lambda x: x/(x**2+1), False, False)
        print "Numerical Result: ", ans, " (required %s)"%k0(1)
        assert np.isclose(ans,k0(1),rtol=1e-3)

    # The following tests follow from the table on Wikipedia
    def test_x2(self):
        """
        Test f(x) = x^2, nu=0
        """
        ht = HankelTransform(nu=0, N=100, h=10**-1.5)

        ans = ht.integrate(lambda x: x**2, False, False)
        print "Numerical Result: ", ans, " (required -1)"
        assert np.isclose(ans,-1,rtol=1e-3)

    def test_x4(self):
        ht = HankelTransform(nu=0, N=150, h=10**-1.5)
        ans = ht.integrate(lambda x: x**4, False, False)
        print "Numerical Result: ", ans, " (required 9)"
        assert np.isclose(ans,9,rtol=1e-3)

    def test_1_on_sqrt_x(self):
        ## NOTE: this is REALLY finnicky!! (check devel/)
        ht = HankelTransform(nu=0, N=160, h=10**-3.5)
        ans = ht.integrate(lambda x: 1./np.sqrt(x), False, False)
        m = -1.5
        anl = 2**(m+1)*gamma(m/2+1)/gamma(-m/2)

        print "Numerical Result: ", ans, " (required %s)"%anl
        assert np.isclose(ans,anl,rtol=1e-3)

    def test_x_on_sqrt_x2_pz2(self):
        # Note that the number required is highly dependent on z .... smaller z is harder.
        ht = HankelTransform(nu=0, N=50, h=10**-1.3)

        z = 1
        ans = ht.integrate(lambda x: x/np.sqrt(x**2+z**2), False, False)
        anl = np.exp(-z)
        print "Numerical Result: ", ans, " (required %s)"%anl
        assert np.isclose(ans,anl,rtol=1e-3)

    def test_gauss(self):
        z = 2
        ht = HankelTransform(nu=0, N=50, h=0.01)

        ans = ht.integrate(lambda x: x*np.exp(-0.5*z**2*x**2), False, False)
        anl = 1./z**2 * np.exp(-0.5/z**2)
        print "Numerical Result: ", ans, " (required %s)"%anl
        assert np.isclose(ans,anl,rtol=1e-3)


class TestAnalyticIntegrals_VaryingOrder(object):
    def powerlaw(self,s,nu,N,h):
        ht = HankelTransform(nu=nu, N=N, h=h)

        ans = ht.integrate(lambda x: x**(s+1), False, False)
        anl = 2**(s+1)*gamma(0.5*(2+nu+s))/gamma(0.5*(nu-s))

        print "Numerical Result: ", ans, " (required %s)"%anl
        assert np.isclose(ans,anl,rtol=1e-3)

    def test_powerlaw(self):
        trials = [[0, 1, 50, 0.05],
                  [0, 2, 50, 0.05],
                  [0.5, 1, 50, 0.05],
                  [-2, 2, 600, 10**-2.6], # This is pretty finnicky
                  [-0.783, 1, 50, 0.05]]

        for s,nu,N,h in trials:
            yield self.powerlaw, s,nu,N,h

    def gamma_mod(self,s,nu,N,h):
        ht = HankelTransform(nu=nu, N=N, h=h)

        ans = ht.integrate(lambda x: x**(nu-2*s + 1)*gammainc_(s,x**2), False, False)
        anl = 0.5 ** (2*s-nu-1) * gammaincc_(1-s+nu,0.25)

        print "Numerical Result: ", ans, " (required %s)"%anl
        assert np.isclose(ans,anl,rtol=1e-3)

    def test_gamma_mod(self):
        trials = [[0, 1, 50, 0.05],
                  [0, 2, 50, 0.05],
                  [0.5, 1, 50, 0.05],
                  [-2, 2, 600, 10**-2.6], # This is pretty finnicky
                  [-0.783, 1, 50, 0.05],
                  [1.0, 0.5, 150, 0.05],
                  [-1., 0.5, 50, 0.05],]

        for s,nu,N,h in trials:
            yield self.powerlaw, s,nu,N,h
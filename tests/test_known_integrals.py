r"""
This module provides some tests of the integrator to known integrals.

Note, these are not the transformations, just the plain integrals,  :math:`\int_0^\infty f(x) J_\nu(x) dx`

There is a corresponding notebook in devel/ that runs each of these functions through a grid of N and h,
showing the pattern of accuracy. This could be useful for finding the correct numbers to choose for these
for other unknown functions.
"""

import pytest

import numpy as np
from scipy.special import gamma, gammainc, gammaincc, k0

from hankel import HankelTransform


def gammainc_(a, x):
    return gamma(a) * gammainc(a, x)


def gammaincc_(a, x):
    return gamma(a) * gammaincc(a, x)


def test_nu0_f_unity():
    """
    Test f(x) = 1, nu=0

    This test is done in the Ogata (2005) paper, section 5"
    """

    ht = HankelTransform(nu=0, N=50, h=0.1)
    ans = ht.integrate(lambda x: 1, False, False)
    print("Numerical Result: ", ans, " (required %s)" % 1)
    assert np.isclose(ans, 1, rtol=1e-3)


def test_nu0_f_x_on_x2():
    """
    Test f(x) = x/(x**2 + 1), nu=0

    This test is done in the Ogata (2005) paper, section 5"
    """
    ht = HankelTransform(nu=0, N=50, h=10 ** -1.5)

    ans = ht.integrate(lambda x: x / (x ** 2 + 1), False, False)
    print("Numerical Result: ", ans, " (required %s)" % k0(1))
    assert np.isclose(ans, k0(1), rtol=1e-3)


def test_nu0_f_x2():
    """
    Test f(x) = x^2, nu=0

    Result on wikipedia
    """
    ht = HankelTransform(nu=0, N=100, h=10 ** -1.5)

    ans = ht.integrate(lambda x: x ** 2, False, False)
    print("Numerical Result: ", ans, " (required -1)")
    assert np.isclose(ans, -1, rtol=1e-3)


def test_nu0_x4():
    """
    Result on wikipedia
    """
    ht = HankelTransform(nu=0, N=150, h=10 ** -1.5)
    ans = ht.integrate(lambda x: x ** 4, False, False)
    print("Numerical Result: ", ans, " (required 9)")
    assert np.isclose(ans, 9, rtol=1e-3)


def test_nu0_1_on_sqrt_x():
    """
    Result on wikipedia
    """
    # NOTE: this is REALLY finnicky!! (check devel/)
    ht = HankelTransform(nu=0, N=160, h=10 ** -3.5)
    ans = ht.integrate(lambda x: 1.0 / np.sqrt(x), False, False)
    m = -1.5
    anl = 2 ** (m + 1) * gamma(m / 2 + 1) / gamma(-m / 2)

    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)


def test_nu0_x_on_sqrt_x2_pz2():
    """
    Result on wikipedia
    """
    # Note that the number required is highly dependent on z .... smaller z is harder.
    ht = HankelTransform(nu=0, N=50, h=10 ** -1.3)

    z = 1
    ans = ht.integrate(lambda x: x / np.sqrt(x ** 2 + z ** 2), False, False)
    anl = np.exp(-z)
    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)


def test_nu0_f_gauss():
    """
    Result on wikipedia
    """
    z = 2
    ht = HankelTransform(nu=0, N=50, h=0.01)

    ans = ht.integrate(lambda x: x * np.exp(-0.5 * z ** 2 * x ** 2), False, False)
    anl = 1.0 / z ** 2 * np.exp(-0.5 / z ** 2)
    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)


@pytest.mark.parametrize(
    "s, nu, N, h",
    [
        [0, 1, 50, 0.05],
        [0, 2, 50, 0.05],
        [0.5, 1, 50, 0.05],
        [-2, 2, 600, 10 ** -2.6],  # This is pretty finnicky
        [-0.783, 1, 50, 0.05],
    ],
)
def test_nu_varying_powerlaw(s, nu, N, h):
    ht = HankelTransform(nu=nu, N=N, h=h)

    ans = ht.integrate(lambda x: x ** (s + 1), False, False)
    anl = 2 ** (s + 1) * gamma(0.5 * (2 + nu + s)) / gamma(0.5 * (nu - s))

    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)


@pytest.mark.parametrize(
    "s, nu, N, h", [[0.5, 1, 50, 0.05], [0.783, 1, 50, 0.05], [1.0, 0.5, 500, 0.01]],
)
def test_nu_varying_gamma_mod(s, nu, N, h):
    ht = HankelTransform(nu=nu, N=N, h=h)

    ans = ht.integrate(
        lambda x: x ** (nu - 2 * s + 1) * gammainc_(s, x ** 2), False, False
    )
    anl = 0.5 ** (2 * s - nu - 1) * gammaincc_(1 - s + nu, 0.25)

    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)

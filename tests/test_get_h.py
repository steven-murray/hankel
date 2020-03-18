r"""
This module provides some tests of the integrator to known integrals.

Note, these are not the transformations, just the plain integrals,  :math:`\int_0^\infty f(x) J_\nu(x) dx`

In this module, we use the get_h function, rather than choosing N and h,
to ensure that this function works.
"""

import pytest

import numpy as np
from scipy.special import gamma, gammainc, gammaincc, k0

from hankel import HankelTransform, SymmetricFourierTransform, get_h


def gammainc_(a, x):
    return gamma(a) * gammainc(a, x)


def gammaincc_(a, x):
    return gamma(a) * gammaincc(a, x)


@pytest.mark.parametrize(
    "f, anl",
    [
        (lambda x: 1, 1),  # Ogata 05
        (lambda x: x / (x ** 2 + 1), k0(1)),  # Ogata 05
        (lambda x: x ** 2, -1),  # wikipedia
        (lambda x: x ** 4, 9),  # wikipedia
        (
            lambda x: 1.0 / np.sqrt(x),  # wikipedia
            2 ** (-0.5) * gamma(-1.5 / 2 + 1) / gamma(1.5 / 2),
        ),
        (lambda x: x / np.sqrt(x ** 2 + 1 ** 2), np.exp(-1)),  # wikipedia
        (
            lambda x: x * np.exp(-0.5 * 2 ** 2 * x ** 2),  # wikipedia
            1.0 / 2 ** 2 * np.exp(-0.5 / 2 ** 2),
        ),
    ],
)
def test_nu0(f, anl):
    ans = get_h(f=f, nu=0, hstart=0.5, atol=0, rtol=1e-3, maxiter=20)[1]
    print("Numerical Result: {ans} (required {anl})".format(ans=ans, anl=anl))
    # we ensure that the answer is within 5e-3 of the true answer,
    # because the algorithm does not ensure that the result is within that
    # tolerance of the true answer, only the previous iteration.
    assert np.isclose(ans, anl, atol=0, rtol=5e-3)


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
    # For this one we test the transform instead
    ans = get_h(f=lambda x: x ** s, nu=nu, K=1, hstart=0.5, atol=1e-3, rtol=1e-3)[1]

    anl = 2 ** (s + 1) * gamma(0.5 * (2 + nu + s)) / gamma(0.5 * (nu - s))
    print("Numerical Result: {ans} (required {anl})".format(ans=ans, anl=anl))
    assert np.isclose(ans, anl, rtol=1e-3, atol=1e-3)


@pytest.mark.parametrize(
    "s, nu, N, h", [[0.5, 1, 50, 0.05], [0.783, 1, 50, 0.05], [1.0, 0.5, 500, 0.01]],
)
def test_nu_varying_gamma_mod(s, nu, N, h):
    ans = get_h(
        f=lambda x: x ** (nu - 2 * s + 1) * gammainc_(s, x ** 2),
        nu=nu,
        hstart=0.5,
        atol=1e-3,
        rtol=1e-3,
    )[1]

    anl = 0.5 ** (2 * s - nu - 1) * gammaincc_(1 - s + nu, 0.25)

    print("Numerical Result: {ans} (required {anl})".format(ans=ans, anl=anl))
    assert np.isclose(ans, anl, rtol=1e-3, atol=1e-3)


def test_too_large_h():
    with pytest.raises(ValueError):
        get_h(f=lambda x: 1, nu=0, hstart=2)


def test_never_converges_high():
    with pytest.raises(Exception):
        get_h(f=lambda x: x ** 6, nu=0, maxiter=15)


def test_never_converges_low():
    with pytest.raises(Exception):
        get_h(f=lambda x: x ** -6, nu=0, maxiter=15)


def test_xrange():
    x1 = HankelTransform.xrange_approx(h=0.5, nu=0)
    x2 = HankelTransform.xrange_approx(h=0.1, nu=0)

    assert x1[0] > x2[0]
    assert x1[1] < x2[1]

    x1 = SymmetricFourierTransform.xrange_approx(h=0.5, ndim=2)
    x2 = SymmetricFourierTransform.xrange_approx(h=0.1, ndim=2)

    assert x1[0] > x2[0]
    assert x1[1] < x2[1]


def test_final_term_amplitude():
    g1 = HankelTransform.final_term_amplitude(f=lambda x: x ** -4, h=0.5)
    g2 = HankelTransform.final_term_amplitude(f=lambda x: x ** -4, h=0.1)
    assert g1 > g2

    g1 = SymmetricFourierTransform.final_term_amplitude(f=lambda x: x ** -4, h=0.5)
    g2 = SymmetricFourierTransform.final_term_amplitude(f=lambda x: x ** -4, h=0.1)
    assert g1 > g2

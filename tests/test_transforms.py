"""
This module provides some tests of the transform, based on analytic transforms provided by Wikipedia.

Note that most of the tests in here are mirrored in test_known_integrals, in which k is set to 1.

There is a corresponding notebook in devel/ that runs each of these functions through a grid of N and h,
showing the pattern of accuracy. This could be useful for finding the correct numbers to choose for these
for other unknown functions.
"""

import pytest

import numpy as np
from scipy.special import gamma, gammainc, gammaincc

from hankel import HankelTransform


def gammainc_(a, x):
    return gamma(a) * gammainc(a, x)


def gammaincc_(a, x):
    return gamma(a) * gammaincc(a, x)


@pytest.mark.parametrize(
    "s, nu, k, N, h",
    [
        [-2, 1, 0.01, 300, 10 ** -3.2],
        [-2, 1, 1, 300, 10 ** -3.2],
        [-2, 1, 10.0, 300, 10 ** -3.2],
        [-2, 2, 0.01, 200, 10 ** -2.0],
        [-2, 2, 1, 200, 10 ** -2.0],
        [-2, 2, 10.0, 200, 10 ** -2.0],
        [-1, 0, 0.01, 50, 0.05],
        [-1, 0, 1, 50, 0.05],
        [-1, 0, 10.0, 50, 0.05],
        [-1, 1, 0.01, 50, 0.05],
        [-1, 1, 1, 50, 0.05],
        [-1, 1, 10.0, 50, 0.05],
        [-1, 2, 0.01, 50, 0.05],
        [-1, 2, 1, 50, 0.05],
        [-1, 2, 10.0, 50, 0.05],
        [1, 0, 0.01, 150, 10 ** -1.5],
        [1, 0, 1, 150, 10 ** -1.5],
        [1, 0, 10.0, 150, 10 ** -1.5],
        [1, 2, 0.01, 50, 0.05],
        [1, 2, 1, 50, 0.05],
        [1, 2, 10.0, 50, 0.05],
        [2, 1, 0.01, 100, 10 ** -1.5],
        [2, 1, 1, 100, 10 ** -1.5],
        [2, 1, 10.0, 100, 10 ** -1.5],
    ],
)
def test_powerlaw(s, nu, k, N, h):
    """
    Test f(r) = 1/r, nu=0
    """

    ht = HankelTransform(nu=nu, N=N, h=h)
    ans = ht.transform(lambda x: x ** s, k, False, False)
    if nu - s <= 0 and (nu - s) % 2 == 0:
        raise Exception("Can't have a negative integer for gamma")

    anl = (
        2 ** (s + 1) * gamma(0.5 * (2 + nu + s)) / k ** (s + 2) / gamma(0.5 * (nu - s))
    )

    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)


@pytest.mark.parametrize(
    "s, nu, k, N, h",
    [
        [-2, 1, 0.01, 300, 10 ** -3.2],
        [-2, 1, 1, 300, 10 ** -3.2],
        [-2, 1, 10.0, 300, 10 ** -3.2],
        [-2, 2, 0.01, 200, 10 ** -2.0],
        [-2, 2, 1, 200, 10 ** -2.0],
        [-2, 2, 10.0, 200, 10 ** -2.0],
        [-1, 0, 0.01, 50, 0.05],
        [-1, 0, 1, 50, 0.05],
        [-1, 0, 10.0, 50, 0.05],
        [-1, 1, 0.01, 50, 0.05],
        [-1, 1, 1, 50, 0.05],
        [-1, 1, 10.0, 50, 0.05],
        [-1, 2, 0.01, 50, 0.05],
        [-1, 2, 1, 50, 0.05],
        [-1, 2, 10.0, 50, 0.05],
        [1, 0, 0.01, 150, 10 ** -1.5],
        [1, 0, 1, 150, 10 ** -1.5],
        [1, 0, 10.0, 150, 10 ** -1.5],
        [1, 2, 0.01, 50, 0.05],
        [1, 2, 1, 50, 0.05],
        [1, 2, 10.0, 50, 0.05],
        [2, 1, 0.01, 100, 10 ** -1.5],
        [2, 1, 1, 100, 10 ** -1.5],
        [2, 1, 10.0, 100, 10 ** -1.5],
    ],
)
def test_alternative(s, nu, k, N, h):
    """Test alternative hankel definition."""
    ht1 = HankelTransform(nu=nu, N=N, h=h)
    ht2 = HankelTransform(nu=nu, N=N, h=h, alt=True)
    ft1 = ht1.transform(lambda r: r ** s, k, False, False)
    ft2 = ht2.transform(lambda r: r ** (s + 0.5), k, False, False) / k ** 0.5
    print("Numerical Results: ", ft1, " and ", ft2)
    assert np.isclose(ft1, ft2, rtol=1e-3)


@pytest.mark.parametrize(
    "nu, alt",
    [
        [-0.5, False],
        [0, False],
        [0.5, False],
        [1, False],
        [1.5, False],
        [-0.5, True],
        [0, True],
        [0.5, True],
        [1, True],
        [1.5, True],
    ],
)
def test_k_zero(nu, alt):
    """Testing k=0."""
    threshold = -0.5 if alt else 0
    ht = HankelTransform(nu=nu, N=50, h=1e-3, alt=alt)
    ans = ht.transform(lambda r: np.exp(-(r ** 2) / 2), 0, False, False)
    print("Numerical Results: ", ans)
    if nu < threshold:
        assert np.isnan(ans)
    elif nu > threshold:
        assert np.isclose(ans, 0, rtol=1e-3)
    else:
        assert np.isclose(ans, 1, rtol=1e-3)

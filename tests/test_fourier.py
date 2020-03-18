"""
This module provides some tests for the Symmetric Fourier Transforms.

Note that most of the tests in here are mirrored in test_known_integrals, in which k is set to 1.

There is a corresponding notebook in devel/ that runs each of these functions through a grid of N and h,
showing the pattern of accuracy. This could be useful for finding the correct numbers to choose for these
for other unknown functions.
"""

import pytest

import numpy as np
from scipy.special import gamma, gammainc, gammaincc

from hankel import SymmetricFourierTransform


def gammainc_(a, x):
    return gamma(a) * gammainc(a, x)


def gammaincc_(a, x):
    return gamma(a) * gammaincc(a, x)


@pytest.mark.parametrize(
    "s, ndim, k, N, h",
    [
        [-0.5, 1, 0.01, 50, 0.05],
        [-0.5, 1, 1, 50, 0.05],
        [-0.5, 1, 10.0, 50, 0.05],
        [-1, 2, 0.01, 50, 0.05],
        [-1, 2, 1, 50, 0.05],
        [-1, 2, 10.0, 50, 0.05],
        [-1, 4, 0.01, 50, 0.05],
        [-1, 4, 1, 50, 0.05],
        [-1, 4, 10.0, 50, 0.05],
        [-1, 6, 0.01, 50, 0.05],
        [-1, 6, 1, 50, 0.05],
        [-1, 6, 10.0, 50, 0.05],
        [1, 1, 0.01, 150, 10 ** -1.5],
        [1, 1, 1, 150, 10 ** -1.5],
        [1, 1, 10.0, 150, 10 ** -1.5],
        [1, 2, 0.01, 150, 10 ** -1.5],
        [1, 2, 1, 150, 10 ** -1.5],
        [1, 2, 10.0, 150, 10 ** -1.5],
        [1, 3, 0.01, 150, 10 ** -1.5],
        [1, 3, 1, 150, 10 ** -1.5],
        [1, 3, 10.0, 150, 10 ** -1.5],
        [-1, 3, 0.01, 150, 10 ** -1.5],
        [-1, 3, 1, 150, 10 ** -1.5],
        [-1, 3, 10.0, 150, 10 ** -1.5],
    ],
)
def test_powerlaw(s, ndim, k, N, h):
    ht = SymmetricFourierTransform(ndim=ndim, N=N, h=h)
    ans = ht.transform(lambda x: x ** s, k, False, False)

    nu = ndim / 2.0 - 1
    s += nu
    if nu - s <= 0 and (nu - s) % 2 == 0:
        raise Exception("Can't have a negative integer for gamma")

    anl = (
        (2 * np.pi) ** (ndim / 2.0)
        * 2 ** (s + 1)
        * gamma(0.5 * (2 + nu + s))
        / k ** (s + 2)
        / gamma(0.5 * (nu - s))
        / k ** nu
    )

    print("Numerical Result: ", ans, " (required %s)" % anl)
    assert np.isclose(ans, anl, rtol=1e-3)

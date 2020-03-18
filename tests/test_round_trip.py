"""
This module provides tests that ensure that the forward then backward transform ends up
with the original.

For simplicity and because of common use (for me at least), we do this with a broken
powerlaw, which turns over at some x0. This *can* ensure convergence of the integral, as
long as the power-law exponents are chosen correctly (some here are commented out
because they don't fit the convergence criteria at high x).

There is a corresponding notebook in devel/ that runs each of these functions through a
grid of N and h, showing the pattern of accuracy. This could be useful for finding the
correct numbers to choose for these for other unknown functions.
"""

import pytest

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

from hankel import SymmetricFourierTransform

r = np.array([0.1, 1, 10.0])


@pytest.mark.parametrize(
    "s, t, x0, ndim, N, h, r",
    [
        [0, 1, 1, 2, 1100, 10 ** -2.5, r],
        [0, 2, 1, 2, 800, 10 ** -2.6, r],
        [0, 3, 1, 2, 1500, 10 ** -3.5, r[:2]],  # Doesn't work for high k!
        [1, 2, 1, 2, 250, 10 ** -2, r],
        [1, 3, 1, 2, 1200, 0.0008, r],
        [2, 3, 1, 2, 1200, 10 ** -2.5, r],
        [0, 1, 10, 2, 1200, 10 ** -2.6, r[1:]],  # Doesn't work for low k
        [0, 2, 10, 2, 250, 10 ** -2.0, r[1:]],  # Doesn't work for low k
        [0, 3, 10, 2, 1000, 10 ** -3, r],
        [1, 2, 10, 2, 1100, 10 ** -2.5, r],
        [1, 3, 10, 2, 250, 10 ** -2.2, r[1:]],  # Doesn't work for low k
        [2, 3, 10, 2, 700, 10 ** -2.3, r[1:]],  # Doesn't work for low k
        [0, 2, 1, 3, 700, 10 ** -2.4, r],
        [0, 3, 1, 3, 250, 10 ** -2, r],
        [1, 3, 1, 3, 2000, 10 ** -2.9, r],
        [0, 2, 10, 3, 2200, 10 ** -2.8, r],
        [0, 3, 10, 3, 400, 10 ** -2.4, r],
        [1, 3, 10, 3, 5000, 10 ** -3, r],  # Bad at low k
    ],
)
def test_roundtrip_broken_pl(s, t, x0, ndim, N, h, r):
    def f(x):
        return x ** s / (x ** t + x0)

    ht = SymmetricFourierTransform(ndim=ndim, N=N, h=h)

    k = np.logspace(
        np.log10(ht.x.min() / r.max() / 10.0),
        np.log10(10 * ht.x.max() / r.min()),
        1000,
    )
    resk = ht.transform(f, k, False)
    spl = Spline(k, resk)
    res = ht.transform(spl, r, False, inverse=True)

    print("Analytic: ", f(r), ", Numerical: ", res)
    assert np.allclose(res, f(r), atol=1e-2)

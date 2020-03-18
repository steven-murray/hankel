# -*- coding: utf-8 -*-
r"""Tools for Hankel transformations."""

import numpy as np
from mpmath import fp as mpm
from scipy.special import gamma, j0, j1, jn
from scipy.special import jn_zeros as _jn_zeros
from scipy.special import jv, yv

SRPI2 = np.sqrt(np.pi / 2.0)


def psi(t):
    """Compute the variable transform from Ogata 2005."""
    return t * np.tanh(np.pi * np.sinh(t) / 2)


def d_psi(t):
    """Compute the derivative of the variable transform from Ogata 2005."""
    t = np.array(t, dtype=float)
    a = np.ones_like(t)
    mask = t < 6
    t = t[mask]
    a[mask] = (np.pi * t * np.cosh(t) + np.sinh(np.pi * np.sinh(t))) / (
        1.0 + np.cosh(np.pi * np.sinh(t))
    )
    return a


def weight(nu, zeros):
    """Get weights for the summation in the hankel transformation."""
    return yv(nu, np.pi * zeros) / kernel(np.pi * zeros, nu + 1)


def roots(N, nu):
    """Get the first N Roots of the Bessel J(nu) functions divided by pi."""
    if np.isclose(nu, np.floor(nu)):
        return _jn_zeros(nu, N) / np.pi
    if np.isclose(nu, 0.5):
        # J(0.5) is just sqrt(2/(x*pi))*sin(x)
        return np.arange(1, N + 1)
    if np.isclose(nu, -0.5):
        # J(-0.5) is just sqrt(2/(x*pi))*cos(x)
        return np.arange(1, N + 1) - 0.5
    return np.array([mpm.besseljzero(nu, i + 1) for i in range(N)]) / np.pi


def j_lim(nu):
    """
    Compute the timit factor of Bessel J(nu, 0) = 0.5 ** nu / Gamma(nu + 1).

    Parameters
    ----------
    nu : float
        Order of the Bessel function.

    Returns
    -------
    float
        The factor.
    """
    return 0.5 ** nu / gamma(nu + 1)


def kernel(x, nu, alt=False):
    """
    Compute kernel functions for the hankel transformation.

    J(nu, x) or for alt=True: J(nu, x) * sqrt(x).

    Parameters
    ----------
    x : array-like
        input values.
    nu : int or float
        order of the bessel function.
    alt : bool, optional
        Whether the alternative defintion of the hankel transform should be
        used: J(nu, x)*sqrt(x). The default is False.

    Returns
    -------
    array-like
        The needed function for the hankel transformation.

    Notes
    -----
    J(nu, x) is approximately (x/2)^nu / Gamma(nu+1) for small x.
    """
    if alt:
        if np.isclose(nu, 0):
            return j0(x) * np.sqrt(x)
        if np.isclose(nu, 1):
            return j1(x) * np.sqrt(x)
        if np.isclose(nu, 0.5):  # J[0.5] = sqrt(2/(x*pi))*sin(x)
            return np.sin(x) / SRPI2
        if np.isclose(nu, -0.5):  # J[-0.5] = sqrt(2/(x*pi))*cos(x)
            return np.cos(x) / SRPI2
        if np.isclose(nu, np.floor(nu)):
            return jn(int(nu), x) * np.sqrt(x)
        return jv(nu, x) * np.sqrt(x)
    if np.isclose(nu, 0):
        return j0(x)
    if np.isclose(nu, 1):
        return j1(x)
    if np.isclose(nu, np.floor(nu)):
        return jn(int(nu), x)
    return jv(nu, x)


def safe_power(x, p):
    """
    Safely calculate x**p.

    Parameters
    ----------
    x : array-like
        value.
    p : float
        exponent.

    Returns
    -------
    array-like
        The result x**p.
    """
    return np.ones_like(x) if np.isclose(p, 0) else np.array(x ** p)


def get_x(h, zeros):
    """Get the arguments for the integrand."""
    return np.pi * psi(h * zeros) / h


def fourier_norm(a, b, ndim, inverse=False):
    r"""Calculate fourier-pair normalisations."""
    if inverse:
        return np.sqrt(np.abs(b) / (2 * np.pi) ** (1 + a)) ** ndim
    return np.sqrt(np.abs(b) / (2 * np.pi) ** (1 - a)) ** ndim


def dim_to_nu(ndim):
    """Calculate nu for the Hankel Transformation."""
    # keep int-type in python 3
    if np.isclose(ndim % 2, 0):
        return int(ndim) // 2 - 1
    return ndim / 2.0 - 1


def get_h(
    f,
    nu,
    K=None,
    cls=None,
    hstart=0.05,
    hdecrement=2,
    atol=1e-3,
    rtol=1e-3,
    maxiter=15,
    inverse=False,
):
    """
    Determine the largest value of h which gives a converged solution.

    Parameters
    ----------
    f : callable
        The function to be integrated/transformed.
    nu : float
        Either the order of the transformation, or the number of dimensions
        (if `cls` is a :class:`SymmetricFourierTransform`)
    K : float or array-like, optional
        The scale(s) of the transformation.
        If None, assumes an integration over f(x)J_nu(x) is desired.
        It is recommended to use a down-sampled K
        for this routine for efficiency. Often a min/max is enough.
    cls : :class:`HankelTransform` subclass, optional
        Either :class:`HankelTransform` or a subclass,
        specifying the type of transformation to do on `f`.
    hstart : float, optional
        The starting value of h.
    hdecrement : float, optional
        How much to divide h by on each iteration.
    atol, rtol : float, optional
        The tolerance parameters, passed to `np.isclose`,
        defining the stopping condition.
    maxiter : int, optional
        Maximum number of iterations to perform.
    inverse : bool, optional
        Whether to treat as an inverse transformation.

    Returns
    -------
    h : float
        The h value at which the solution converges.
    res : scalar or tuple
        The value of the integral/transformation using the returned h --
        if a transformation, returns results at K.
    N : int
        The number of nodes necessary in the final calculation.
        While each iteration uses N=pi/h, the returned N checks
        whether nodes are numerically zero above some threshold.

    Notes
    -----
    This function is not completely general. The function `f` is assumed to be
    reasonably smooth and non-oscillatory.

    The idea is to use successively smaller values of *h*, with N=pi/h on each
    iteration, until the result betweeniterations becomes stable.
    """
    # prevent circular imports
    from hankel.hankel import HankelTransform

    cls = HankelTransform if cls is None else cls

    if hstart >= 1:
        raise ValueError("h should never be greater than unity")

    # First, ensure that *some* of the values are non-zero
    i = 0
    while (
        np.any(
            np.all(
                cls(nu, h=hstart, N=int(np.pi / hstart))._get_series(
                    f, 1 if K is None else K
                )
                == 0,
                axis=-1,
            )
        )
        and i < maxiter
    ):
        hstart /= hdecrement
        i += 1

    if i == maxiter:
        raise Exception("Maxiter reached while checking for non-zero values")

    if K is None:  # Do a normal integral of f(x)J_nu(x)
        K = 1

        def getres(h):
            """Get the result."""
            return cls(nu, h=h, N=int(np.pi / h)).transform(
                lambda x: f(x) / x, k=K, ret_err=False, inverse=inverse
            )

    else:  # Do a transform at k=K

        def getres(h):
            """Get the result."""
            return cls(nu, h=h, N=int(np.pi / h)).transform(
                f, k=K, ret_err=False, inverse=inverse
            )

    res = getres(hstart)
    res2 = 2 * res + 10

    while not np.allclose(res, res2, atol=atol, rtol=rtol) and i < maxiter:
        i += 1
        hstart /= hdecrement
        res2 = 1 * res
        res = getres(hstart)

    if i == maxiter:
        raise Exception("Maxiter reached while checking convergence")

    # Can do some more trimming of N potentially, by seeing where f(x)~0.
    def consecutive(data, stepsize=1):
        """Split into arrays of consecutive zeros."""
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)

    hstart *= hdecrement

    x = cls(nu, h=hstart, N=int(np.pi / hstart)).x
    lastk = np.where(f(x / np.max(K)) == 0)[0]
    if len(lastk) > 1:
        # if there are any that are zero,
        # and if there are more than 1 in a row
        # (otherwise might just be oscillatory)
        lastk = consecutive(lastk)  # split into arrays of consecutive zeros
        if len(lastk[-1]) == 1:
            lastk = int(np.pi / hstart)
        else:
            lastk = lastk[-1][0]
    else:  # otherwise set back to N
        lastk = int(np.pi / hstart)

    return hstart, res2, lastk

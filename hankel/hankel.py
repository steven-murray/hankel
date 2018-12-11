r"""
General quadrature method for Hankel transformations.

Based on the algorithm provided in
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
Publications of the Research Institute for Mathematical Sciences,
vol. 41, no. 4, pp. 949-970, 2005.
"""

# TODO: Suppress warnings on overflows
# TODO: Write tests
# TODO: Profile.
from __future__ import division

import numpy as np
from mpmath import fp as mpm
from scipy.special import j0, j1, jn_zeros as _jn_zeros, jn, yv, jv


class HankelTransform(object):
    r"""
    The basis of the Hankel Transformation algorithm by Ogata 2005.

    This algorithm is used to solve the equation
    :math:`\int_0^\infty f(x) J_\nu(x) dx`
    where :math:`J_\nu(x)` is a Bessel function of the first kind of order
    :math:`nu`, and :math:`f(x)` is an arbitrary (slowly-decaying) function.

    The algorithm is presented in
    H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
    Publications of the Research Institute for Mathematical Sciences,
    vol. 41, no. 4, pp. 949-970, 2005.

    This class provides a method for directly performing this integration,
    and also for doing a Hankel Transform.

    Parameters
    ----------

    nu : int or 0.5, optional, default = 0
        The order of the bessel function (of the first kind) J_nu(x)

    N : int, optional, default = 3.2/`h`
        The number of nodes in the calculation. Generally this must increase
        for a smaller value of the step-size h. Default value is based on where
        the series will truncate according
        to the double-exponential convergence to
        the roots of the Bessel function.

    h : float, optional, default = 0.1
        The step-size of the integration.
    """

    def __init__(self, nu=0, N=None, h=0.05):

        if N is None:
            N = int(3.2 / h)

        if not np.isscalar(N):
            raise ValueError("N must be a scalar")
        if not np.isscalar(h):
            raise ValueError("h must be a scalar")
        if not np.isscalar(nu):
            raise ValueError("nu must be a scalar")

        self._nu = nu
        self._h = h
        self._zeros = self._roots(N)
        self.x = self._x(h)
        self.j = self._j(self.x)
        self.w = self._weight()
        self.dpsi = self._d_psi(h * self._zeros)

        # Some quantities only useful in the FourierTransform
        self._x_power = 1
        self._k_power = 2

    def _psi(self, t):
        y = np.sinh(t)
        return t * np.tanh(np.pi * y / 2)

    def _d_psi(self, t):
        a = (np.pi * t * np.cosh(t) + np.sinh(np.pi * np.sinh(t))) / (
            1.0 + np.cosh(np.pi * np.sinh(t))
        )
        a[np.isnan(a)] = 1.0
        return a

    def _weight(self):
        return yv(self._nu, np.pi * self._zeros) / self._j1(
            np.pi * self._zeros
        )

    def _roots(self, N):
        if isinstance(self._nu, int):
            return _jn_zeros(self._nu, N) / np.pi
        elif np.isclose(self._nu, 0.5):
            # J[0.5] = sqrt(2/(x*pi))*sin(x)
            return np.arange(1, N + 1)
        elif np.isclose(self._nu, -0.5):
            # J[-0.5] = sqrt(2/(x*pi))*cos(x)
            return np.arange(1, N + 1) - 0.5
        else:
            return (
                np.array([mpm.besseljzero(self._nu, i + 1) for i in range(N)])
                / np.pi
            )

    def _j(self, x):
        if self._nu == 0:
            return j0(x)
        elif self._nu == 1:
            return j1(x)
        elif isinstance(self._nu, int):
            return jn(self._nu, x)
        else:
            return jv(self._nu, x)

    def _j1(self, x):
        if self._nu == -1:
            return j0(x)
        elif self._nu == 0:
            return j1(x)
        elif isinstance(self._nu, int):
            return jn(self._nu + 1, x)
        else:
            return jv(self._nu + 1, x)

    def _x(self, h):
        return np.pi * self._psi(h * self._zeros) / h

    def _f(self, f, x):
        return f(x)

    @staticmethod
    def _k(k):
        return k

    @staticmethod
    def _norm(self, inverse=False):
        r"""
        Scalar normalisation of the transform. Identically 1.
        """
        return 1

    def _get_series(self, f, k=1):
        fres = (
            self._f(f, np.divide.outer(self.x, k).T) * self.x ** self._x_power
        )
        return np.pi * self.w * fres * self.j * self.dpsi

    def transform(self, f, k=1, ret_err=True, ret_cumsum=False, inverse=False):
        r"""
        Do the Hankel-transform of the function f.

        Parameters
        ----------
        f : callable
            A function of one variable, representing :math:`f(x)`

        ret_err : boolean, optional, default = True
            Whether to return the estimated error

        ret_cumsum : boolean, optional, default = False
            Whether to return the cumulative sum

        Returns
        -------
        ret : array-like
            The Hankel-transform of f(x) at the provided k. If
            `k` is scalar, then this will be scalar.

        err : array-like
            The estimated error of the approximate integral, at every `k`.
            It is merely the last term in the sum.
            Only returned if `ret_err=True`.

        cumsum : array-like
            The total cumulative sum, for which the last term
            is itself the transform.
            One can use this to check whether the integral is converging.
            Only returned if `ret_cumsum=True`


        Notes
        -----
        The Hankel transform is defined as

        .. math:: F(k) = \int_0^\infty r f(r) J_\nu(kr) dr.

        The inverse transform is identical (swapping *k* and *r* of course).
        """
        # The following allows for a re-scaling of k when doing FT's.
        k = self._k(k)

        # The following is the scalar normalisation of the transform
        # The basic transform has a norm of 1.
        # But when doing FT's, this depends on the dimensionality.
        norm = self._norm(inverse)

        # The following renormalises by the fourier dual to some power
        knorm = k ** self._k_power

        summation = self._get_series(f, k)
        ret = norm * np.sum(summation, axis=-1) / knorm

        if ret_err:
            err = norm * np.take(summation, -1, axis=-1) / knorm
        if ret_cumsum:
            cumsum = norm * np.cumsum(summation, axis=-1).T / knorm

        if ret_err and ret_cumsum:
            return ret, err, cumsum
        elif ret_err:
            return ret, err
        elif ret_cumsum:
            return ret, cumsum
        else:
            return ret

    def integrate(self, f, ret_err=True, ret_cumsum=False):
        r"""
        Do the Hankel-type integral of the function f.

        This is *not* the Hankel transform, but rather the simplified
        integral, :math:`\int_0^\infty f(x) J_\nu(x) dx` , equivalent to the
        transform of :math:`f(r)/r` at *k=1*.

        Parameters
        ----------
        f : callable
            A function of one variable, representing :math:`f(x)`

        ret_err : boolean, optional, default = True
            Whether to return the estimated error

        ret_cumsum : boolean, optional, default = False
            Whether to return the cumulative sum
        """
        return self.transform(
            f=lambda x: f(x) / x,
            k=1,
            ret_err=ret_err,
            ret_cumsum=ret_cumsum,
            inverse=False,
        )

    def xrange(self, k=1):
        """
        Tuple giving (min,max) x value evaluated by f(x).

        Parameters
        ----------
        k : array-like, optional
            Scales for the transformation. Leave as 1 for an integral.

        See Also
        --------
        See :meth:`xrange_approx` for an approximate version of this method
        which is a classmethod.
        """
        return np.array([self.x.min() / np.max(k), self.x.max() / np.min(k)])

    @classmethod
    def xrange_approx(cls, h, nu, k=1):
        """
        Tuple giving approximate (min,max) x value evaluated by f(x/k).

        Operates under the assumption that N = 3.2/h.

        Parameters
        ----------
        h : float
            The resolution parameter of the Hankel integration
        nu : float
            Order of the integration/transform
        k : array-like, optional
            Scales for the transformation. Leave as 1 for an integral.

        See Also
        --------
        See :meth:`xrange` (instance method) for the actual x-range
        under a given choice of parameters.
        """
        r = mpm.besseljzero(nu, 1) / np.pi
        return np.array([np.pi ** 2 * h * r ** 2 / 2 / k, np.pi * 3.2 / h / k])

    @classmethod
    def G(cls, f, h, k=None, *args, **kwargs):
        """
        The absolute value of the non-oscillatory
        of the summed series' last term, up to a scaling constant.

        This can be used to get the sign of the slope of G with h.

        Parameters
        ----------
        f : callable
            The function to integrate/transform
        h : float
            The resolution parameter of the hankel integration
        k : float or array-like, optional
            The scale at which to evaluate the transform.
            If None, assume an integral.

        Returns
        -------
        The value of G.
        """
        if k is None:
            return np.sqrt(2 * h / 3.2) * f(3.2 * np.pi / h)
        else:
            return np.sqrt(3.2 / (2 * h)) * f(3.2 * np.pi / h / k)

    @classmethod
    def deltaG(cls, f, h, *args, **kwargs):
        "The slope (up to a constant) of the last term of the series with h"
        return cls.G(f, h, *args, **kwargs) - cls.G(
            f, h / 1.1, *args, **kwargs
        )


class SymmetricFourierTransform(HankelTransform):
    r"""
    Determine the Fourier Transform of a radially symmetric function
    in arbitrary dimensions.

    Parameters
    ----------
    ndim : int
        Number of dimensions the transform is in.

    a, b : float, default 1
        This pair of values defines the Fourier convention used
        (see Notes below for details)

    N : int, optional
        The number of nodes in the calculation. Generally this must increase
        for a smaller value of the step-size h.

    h : float, optional
        The step-size of the integration.

    Notes
    -----
    We allow for arbitrary Fourier convention, according to the scheme in
    http://mathworld.wolfram.com/FourierTransform.html.
    That is, we define the forward and inverse *n*-dimensional transforms
    respectively as

    .. math:: F(k) = \sqrt{\frac{|b|}{(2\pi)^{1-a}}}^n
              \int f(r) e^{i b\mathbf{k}\cdot\mathbf{r}} d^n\mathbf{r}

    and

    .. math:: f(r) = \sqrt{\frac{|b|}{(2\pi)^{1+a}}}^n
              \int F(k) e^{-i b\mathbf{k}\cdot\mathbf{r}} d^n \mathbf{k}.

    By default, we set both *a* and *b* to 1,
    so that the forward transform has a normalisation of unity.

    In this general sense, the forward and inverse Hankel transforms are
    respectively

    .. math:: F(k) = \sqrt{\frac{|b|}{(2\pi)^{1-a}}}^n
              \frac{(2\pi)^{n/2}}{(bk)^{n/2-1}}
              \int_0^\infty r^{n/2-1} f(r) J_{n/2-1}(bkr) r dr

    and

    .. math:: f(r) = \sqrt{\frac{|b|}{(2\pi)^{1+a}}}^n
              \frac{(2\pi)^{n/2}}{(br)^{n/2-1}}
              \int_0^\infty k^{n/2-1} f(k) J_{n/2-1}(bkr) k dk.

    """

    def __init__(self, ndim=2, a=1, b=1, N=200, h=0.05):
        # keep int-type in python 3
        if np.isclose(ndim % 2, 0):
            nu = int(ndim) // 2 - 1
        else:
            nu = ndim / 2.0 - 1

        self.ndim = ndim
        self.fourier_norm_a = a
        self.fourier_norm_b = b

        super(SymmetricFourierTransform, self).__init__(nu=nu, N=N, h=h)

        self._x_power = self.ndim / 2.0
        self._k_power = self.ndim

    def _fourier_norm(self, inverse=False):
        r"""
        Calculate fourier-pair normalisations.

        See class documentation for details.
        """
        if inverse:
            return (
                np.sqrt(
                    np.abs(self.fourier_norm_b)
                    / (2 * np.pi) ** (1 + self.fourier_norm_a)
                )
                ** self.ndim
            )
        else:
            return (
                np.sqrt(
                    np.abs(self.fourier_norm_b)
                    / (2 * np.pi) ** (1 - self.fourier_norm_a)
                )
                ** self.ndim
            )

    def _norm(self, inverse=False):
        r"""
        The scalar normalisation of the transform,
        taking into account Fourier conventions and a possible inversion.
        """
        return (2 * np.pi) ** (self.ndim / 2.0) * self._fourier_norm(inverse)

    def transform(self, f, k, *args, **kwargs):
        r"""
        Do the *n*-symmetric Fourier transform of the function f.

        Parameters and returns are precisely the same as
        :meth:`HankelTransform.transform`.

        Notes
        -----
        The *n*-symmetric fourier transform is defined
        in terms of the Hankel transform as

        .. math:: F(k) = \frac{(2\pi)^{n/2}}{k^{n/2-1}}
                  \int_0^\infty r^{n/2-1} f(r) J_{n/2-1}(kr)r dr.

        The inverse transform has an inverse normalisation.
        """
        k = self.fourier_norm_b * k
        return super(SymmetricFourierTransform, self).transform(
            f, k, *args, **kwargs
        )

    @classmethod
    def xrange_approx(cls, h, ndim, k=1):
        """
        Tuple giving approximate (min,max) x value evaluated by f(x/k).

        Operates under the assumption that N = 3.2/h.

        Parameters
        ----------
        h : float
            The resolution parameter of the Hankel integration

        ndim : float
            Number of dimensions of the transform.

        k : array-like, optional
            Scales for the transformation. Leave as 1 for an integral.

        See Also
        --------
        See :meth:`xrange` (instance method)
        for the actual x-range under a given choice of parameters.
        """
        return HankelTransform.xrange_approx(h, ndim / 2.0 - 1, k)

    @classmethod
    def G(self, f, h, k=None, ndim=2):
        """
        The absolute value of the non-oscillatory
        of the summed series' last term, up to a scaling constant.

        This can be used to get the sign of the slope of G with h.

        Parameters
        ----------
        f : callable
            The function to integrate/transform
        h : float
            The resolution parameter of the hankel integration
        k : float or array-like, optional
            The scale at which to evaluate the transform.
            If None, assume an integral.
        ndim : float
            The number of dimensions of the transform

        Returns
        -------
        The value of G.
        """
        if k is None:
            return HankelTransform.G(f, h, k)
        else:
            return (3.2 / h) ** ((ndim - 1) / 2.0) * f(3.2 * np.pi / h / k)


def get_h(
    f,
    nu,
    K=None,
    cls=HankelTransform,
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
        While each iteration uses N=3.2/h, the returned N checks
        whether nodes are numerically zero above some threshold.

    Notes
    -----
    This function is not completely general.
    The function `f` is assumed to be reasonably smooth and non-oscillatory.

    """

    # First, ensure that *some* of the values are non-zero
    i = 0
    while (
        np.any(
            np.all(
                cls(nu, h=hstart, N=int(3.2 / hstart))._get_series(
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

    if K is None:  # Do a normal integral of f(x)J_nu(x)

        # First ensure that the derivative of G(h) is negative
        while cls.deltaG(f, hstart) > 0:
            hstart /= hdecrement

        K = 1

        def getres(h):
            """dummy function to get the result"""
            return cls(nu, h=h, N=int(3.2 / h)).transform(
                lambda x: f(x) / x, k=K, ret_err=False, inverse=inverse
            )

    else:  # Do a transform at k=K
        while np.any(cls.deltaG(f, hstart, K, nu) > 0):
            hstart /= hdecrement

        def getres(h):
            """dummy function to get the result"""
            return cls(nu, h=h, N=int(3.2 / h)).transform(
                f, k=K, ret_err=False, inverse=inverse
            )

    res = getres(hstart)
    res2 = 2 * res + 10
    i = 0

    while not np.allclose(res, res2) and i < maxiter:
        i += 1
        hstart /= hdecrement
        res2 = 1 * res
        res = getres(hstart)

    if i == maxiter:
        raise Exception("Maxiter reached")

    # Can do some more trimming of N potentially, by seeing where f(x)~0.
    def consecutive(data, stepsize=1):
        """split into arrays of consecutive zeros"""
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)

    hstart *= hdecrement

    x = cls(nu, h=hstart, N=int(3.2 / hstart)).x
    lastk = np.where(f(x / np.max(K)) == 0)[0]
    if len(lastk) > 1:
        # if there are any that are zero,
        # and if there are more than 1 in a row
        # (otherwise might just be oscillatory)
        lastk = consecutive(lastk)  # split into arrays of consecutive zeros
        if len(lastk[-1]) == 1:
            lastk = int(3.2 / hstart)
        else:
            lastk = lastk[-1][0]
    else:  # otherwise set back to N
        lastk = int(3.2 / hstart)

    return hstart, res2, lastk

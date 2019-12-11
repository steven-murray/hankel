r"""
General quadrature method for Hankel transformations.

Based on the algorithm provided in
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
Publications of the Research Institute for Mathematical Sciences,
vol. 41, no. 4, pp. 949-970, 2005.
"""

from __future__ import division, absolute_import
from builtins import super

import numpy as np
from mpmath import fp as mpm
from scipy.integrate import quad
from hankel.tools import (
    d_psi,
    kernel,
    weight,
    roots,
    get_x,
    fourier_norm,
    dim_to_nu,
    safe_power,
    j_lim,
)


class HankelTransform(object):
    r"""
    The basis of the Hankel Transformation algorithm by Ogata 2005.

    This algorithm is used to solve the equation
    :math:`\int_0^\infty f(x) J_\nu(x) dx`
    where :math:`J_\nu(x)` is a Bessel function of the first kind of order
    :math:`\nu`, and :math:`f(x)` is an arbitrary (slowly-decaying) function.

    The algorithm is presented in
    H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
    Publications of the Research Institute for Mathematical Sciences,
    vol. 41, no. 4, pp. 949-970, 2005.

    This class provides a method for directly performing this integration,
    and also for doing a Hankel Transform.

    Parameters
    ----------
    nu : scalar, optional
        The order of the bessel function (of the first kind) J_nu(x)
    N : int, optional, default = ``pi/h``
        The number of nodes in the calculation. Generally this must increase
        for a smaller value of the step-size h. Default value is based on where
        the series will truncate according
        to the double-exponential convergence to
        the roots of the Bessel function.
    h : float, optional
        The step-size of the integration.
    alt : bool, optional
        Whether to use the alternative definition of the hankel transform.
        should be used. Default: False
    """

    def __init__(self, nu=0, N=None, h=0.05, alt=False):

        N = int(np.pi / h) if N is None else N
        if not np.isscalar(N):
            raise ValueError("N must be a scalar")
        if not np.isscalar(h):
            raise ValueError("h must be a scalar")
        if not np.isscalar(nu):
            raise ValueError("nu must be a scalar")
        if nu < -0.5:
            raise ValueError("nu must be at least -1/2")

        self._nu = nu
        self._h = h
        self._zeros = roots(N, nu)
        self.x = get_x(h, self._zeros)
        self.kernel = kernel(self.x, nu, alt)
        self.w = weight(nu, self._zeros)
        self.dpsi = d_psi(h * self._zeros)
        self.alt = alt

        # Some quantities only useful in the FourierTransform
        self._r_power = 0 if alt else 1
        self._k_power = 0

    @property
    def nu(self):
        """Order of the hankel transform."""
        return self._nu

    def _f(self, f, r):
        """Correct integrand depending on transformation definition."""
        return f(r) * safe_power(r, self._r_power)

    def _k(self, k):
        return np.array(k)

    def _norm(self, inverse=False):
        r"""Scalar normalisation of the transform. Identically 1."""
        return 1

    def _get_series(self, f, k=1):
        with np.errstate(divide="ignore"):  # numpy safely divids by 0
            args = np.divide.outer(self.x, k).T  # x = r*k
        fres = self._f(f, args)
        return np.pi * self.w * fres * self.kernel * self.dpsi

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

        Or in the alternative case with ``alt=True``:

        .. math:: F(k) = \int_0^\infty f(r) \sqrt{kr} J_\nu(kr) dr.

        The inverse transform is identical (swapping *k* and *r* of course).
        """
        # The following allows for a re-scaling of k when doing FT's.
        k_scalar = np.isscalar(k)
        k = self._k(k)
        # k = zero here
        k_0 = np.isclose(k, 0)

        # The following is the scalar normalisation of the transform
        # The basic transform has a norm of 1.
        # But when doing FT's, this depends on the dimensionality.
        norm = self._norm(inverse)

        # The following renormalises by the fourier dual to some power
        knorm = safe_power(k, self._k_power + 1)

        # set replace k=0 with 1 in the following
        knorm[k_0] = 1
        k[k_0] = 1

        summation = self._get_series(f, k)
        ret = np.array(norm * np.sum(summation, axis=-1) / knorm)

        # care about k = 0
        if np.any(k_0):
            # limit of J(nu, 0) considering powers of k
            alt_pow = 0.5 if self.alt else 0  # in alt. def sqrt(rk) involved
            lim_r_pow = self._r_power + alt_pow + self.nu
            nu_th = self._k_power - 1 - alt_pow  # threshold
            lim = j_lim(self.nu)

            def integrand(r):
                return f(r) * safe_power(r, lim_r_pow)

            if np.isclose(self.nu, nu_th):
                ret[k_0] = quad(integrand, 0, np.inf)[0] * lim * norm
            elif self.nu < nu_th:
                ret[k_0] = np.nan
            else:
                ret[k_0] = 0

        if k_scalar:
            ret = ret.item()

        if ret_err:
            err = norm * np.take(summation, -1, axis=-1) / knorm
        if ret_cumsum:
            cumsum = norm * np.cumsum(summation, axis=-1).T / knorm

        if ret_err and ret_cumsum:
            return ret, err, cumsum
        if ret_err:
            return ret, err
        if ret_cumsum:
            return ret, cumsum
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
        if self.alt:
            func = lambda x: f(x) / np.sqrt(x)
        else:
            func = lambda x: f(x) / x
        return self.transform(
            f=func,
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
        xrange: the actual x-range under a given choice of parameters.
        """
        r = mpm.besseljzero(nu, 1) / np.pi
        return np.array(
            [np.pi ** 2 * h * r ** 2 / 2 / k, np.pi * np.pi / h / k]
        )

    @classmethod
    def G(cls, f, h, k=None, *args, **kwargs):
        """
        Info about the last term in the series.

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
            return np.sqrt(2 * h / np.pi) * f(np.pi * np.pi / h)
        return np.sqrt(np.pi / (2 * h)) * f(np.pi * np.pi / h / k)

    @classmethod
    def deltaG(cls, f, h, *args, **kwargs):
        """Slope (up to a constant) of the last term of the series with h."""
        return cls.G(f, h, *args, **kwargs) - cls.G(
            f, h / 1.1, *args, **kwargs
        )


class SymmetricFourierTransform(HankelTransform):
    r"""
    Fourier Transform of a radially symmetric function in arbitrary dimensions.

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
    alt : bool, optional
        State if the alternative definition of the hankel transform
        should be used. Default: False

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

    def __init__(self, ndim=2, a=1, b=1, N=None, h=0.05, alt=True):

        super().__init__(nu=dim_to_nu(ndim), N=N, h=h, alt=alt)

        self.ndim = ndim
        self.fourier_norm_a = a
        self.fourier_norm_b = b
        self._r_power = (ndim - 1) / 2.0 if alt else ndim / 2.0
        self._k_power = (ndim - 1) / 2.0 if alt else ndim / 2.0 - 1

    def _norm(self, inverse=False):
        r"""
        Scalar normalisation of the transform.

        Taking into account Fourier conventions and a possible inversion.
        """
        return (2 * np.pi) ** (self.ndim / 2.0) * fourier_norm(
            self.fourier_norm_a, self.fourier_norm_b, self.ndim, inverse
        )

    def _k(self, k):
        """Substitution for k."""
        return np.array(self.fourier_norm_b * np.array(k))

    @classmethod
    def xrange_approx(cls, h, ndim, k=1):
        """
        Tuple giving approximate (min,max) x value evaluated by f(x/k).

        Operates under the assumption that N = pi/h.

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
        xrange:  the actual x-range under a given choice of parameters.
        """
        return HankelTransform.xrange_approx(h, dim_to_nu(ndim), k)

    @classmethod
    def G(cls, f, h, k=None, ndim=2):
        """
        Info about the last term in the series.

        The absolute value of the non-oscillatory part
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
        fmax = f(cls.xrange_approx(h, ndim, k)[-1])
        return (np.pi / h) ** ((ndim - 1) / 2.0) * fmax

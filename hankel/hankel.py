'''
General quadrature method for Hankel transformations.

Based on the algorithm provided in
H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
Publications of the Research Institute for Mathematical Sciences,
vol. 41, no. 4, pp. 949-970, 2005.
'''

# TODO: Suppress warnings on overflows
# TODO: Write tests
# TODO: Profile.

import numpy as np
from mpmath import fp as mpm
from scipy.special import j0, j1, jn_zeros, jn, yv, jv


class HankelTransform(object):
    """
    The basis of the Hankel Transformation algorithm by Ogata 2005.

    This algorithm is used to solve the equation :math:`\int_0^\infty f(x) J_\nu(x) dx`
    where :math:`J_\nu(x)` is a Bessel function of the first kind of order
    :math:`nu`, and :math:`f(x)` is an arbitrary (slowly-decaying) function.

    The algorithm is presented in
    H. Ogata, A Numerical Integration Formula Based on the Bessel Functions,
    Publications of the Research Institute for Mathematical Sciences, vol. 41, no. 4, pp. 949-970, 2005.

    Parameters
    ----------
    nu : int or 0.5, optional, default = 0
        The order of the bessel function (of the first kind) J_nu(x)

    N : int, optional, default = 100
        The number of nodes in the calculation. Generally this must increase
        for a smaller value of the step-size h.

    h : float, optional, default = 0.1
        The step-size of the integration.
    """

    def __init__(self, nu=0, N=200, h=0.05):
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
        self.dpsi = self._d_psi(h*self._zeros)

        # Some quantities only useful in the FourierTransform
        self._x_power = 1
        self._norm = 1
        self._k_power = 2

    def _psi(self, t):
        y = np.sinh(t)
        return t*np.tanh(np.pi*y/2)

    def _d_psi(self, t):
        a = (np.pi*t*np.cosh(t) + np.sinh(np.pi*np.sinh(t)))/(1.0 + np.cosh(np.pi*np.sinh(t)))
        a[np.isnan(a)] = 1.0
        return a

    def _weight(self):
        return yv(self._nu, np.pi*self._zeros)/self._j1(np.pi*self._zeros)

    def _roots(self, N):
        if isinstance(self._nu, int):
            return jn_zeros(self._nu, N)/np.pi
        elif self._nu == 0.5:
            return np.arange(1, N + 1)
        else:
            return np.array([mpm.besseljzero(self._nu, i + 1) for i in range(N)])/np.pi

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
        return np.pi*self._psi(h*self._zeros)/h

    def _f(self, f, x):
        return f(x)

    def transform(self, f, k=1, ret_err=True, ret_cumsum=False, inverse=False):
        """
        Do the Hankel-transform of the function f.

        Parameters
        ----------
        f : callable
            A function of one variable, representing :math:`f(x)`

        ret_err : boolean, optional, default = True
            Whether to return the estimated error

        ret_cumsum : boolean, optional, default = False
            Whether to return the cumulative sum
        """
        if inverse:
            norm = 1./self._norm
        else:
            norm = self._norm

        fres = self._f(f, np.divide.outer(self.x, k).T)*self.x**self._x_power
        summation = np.pi*self.w*fres*self.j*self.dpsi
        ret = norm * np.sum(summation, axis=-1)/k ** self._k_power

        if ret_err:
            err = norm * np.take(summation, -1, axis=-1)/k ** self._k_power
        if ret_cumsum:
            cumsum = norm * np.divide.outer(np.cumsum(summation, axis=-1), k ** self._k_power)

        if ret_err and ret_cumsum:
            return ret, err, cumsum
        elif ret_err:
            return ret, err
        elif ret_cumsum:
            return ret, cumsum
        else:
            return ret

    def integrate(self, f, ret_err=True, ret_cumsum=False):
        """
        Do the Hankel-type integral of the function f.

        This is *not* the Hankel transform, but rather the simplified
        integral, :math:`\int_0^\infty f(x) J_\nu(x) dx`, equivalent to the
        transform of :math:`f(r)/r` at k=1.

        Parameters
        ----------
        f : callable
            A function of one variable, representing :math:`f(x)`

        ret_err : boolean, optional, default = True
            Whether to return the estimated error

        ret_cumsum : boolean, optional, default = False
            Whether to return the cumulative sum
        """
        return self.transform(f=lambda x: f(x)/x, k=1, ret_err=ret_err, ret_cumsum=ret_cumsum, inverse=False)


class SymmetricFourierTransform(HankelTransform):
    """
    Determine the Fourier Transform of a radially symmetric function in arbitrary dimensions.

    See :class:`HankelTransform` for descriptions of inputs. Additional parameters noted below.

    Parameters
    ----------
    ndim : int
        Number of dimensions the transform is in.
    """

    def __init__(self, ndim=2, N=200, h=0.05):
        if ndim%2 == 0:
            nu = ndim/2 - 1
        else:
            nu = ndim/2. - 1

        self.ndim = ndim

        super(SymmetricFourierTransform, self).__init__(nu=nu, N=N, h=h)

        self._x_power = self.ndim/2.
        self._k_power = self.ndim
        self._norm = (2*np.pi) ** (self.ndim/2.)

    # def transform(self, f, k=1, ret_err=True, ret_cumsum=False, inverse=False):
    #     """
    #     Do the fourier transform of f, where f is radially symmetric.
    #
    #     Parameters
    #     ----------
    #     f : callable
    #         A function of one variable, representing :math:`f(x)`
    #
    #     ret_err : boolean, optional, default = True
    #         Whether to return the estimated error
    #
    #     ret_cumsum : boolean, optional, default = False
    #         Whether to return the cumulative sum
    #     """
    #     fres = self._f(f, np.divide.outer(self.x, k).T)*self.x ** (self.ndim/2.)
    #
    #     summation = np.pi*self.w*fres*self.j*self.dpsi
    #
    #     norm = (2*np.pi) ** (self.ndim/2.)
    #     if inverse:
    #         norm = 1./norm
    #
    #     ret = norm*np.sum(summation, axis=-1)/k ** self.ndim
    #
    #     if ret_err:
    #         err = norm*np.take(summation, -1, axis=-1)/k ** self.ndim
    #     if ret_cumsum:
    #         cumsum = norm*np.divide.outer(np.cumsum(summation, axis=-1), k ** self.ndim)
    #
    #     if ret_err and ret_cumsum:
    #         return ret, err, cumsum
    #     elif ret_err:
    #         return ret, err
    #     elif ret_cumsum:
    #         return ret, cumsum
    #     else:
    #         return ret

# class SphericalHankelTransform(HankelTransform):
#     """
#     Perform spherical hankel transforms.
#
#     Defined as :math:`\int_0^\infty f(x) j_\nu(x) dx
#
#     .. Note :: Only does 0th-order transforms currently.
#     """
#     def __init__(self, nu=0, *args, **kwargs):
#         nu += 0.5
#         super(SphericalHankelTransform, self).__init__(nu, *args, **kwargs)
#
#     def _f(self, f, x):
#         return np.sqrt(np.pi / (2 * x)) * f(x)
#
#     def _roots(self, N):
#         if self._nu == 0.5:
#             return (np.arange(N) + 1)
#         else:
#             return super(SphericalHankelTransform, self)._roots(N)
#

'''
This module provides some tests of the API of hankel.

That is, it tests stuff like whether one can input scalars/arrays, and whether the outputs of functions have the correct
number of entries etc.
'''

import inspect
import numpy as np
import os

LOCATION = "/".join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))).split("/")[:-1])
import sys

sys.path.insert(0, LOCATION)

from nose.tools import raises

from hankel import HankelTransform, SymmetricFourierTransform


@raises(ValueError)
def test_array_N():
    ht = HankelTransform(N=np.linspace(100, 200, 10))


@raises(TypeError)
def test_array_h():
    ht = HankelTransform(h=np.logspace(-3, 0, 10))


@raises(ValueError)
def test_array_nu():
    ht = HankelTransform(nu=np.linspace(0, 1, 2))


def test_k_array():
    k = np.logspace(-3, 3, 10)
    ht = HankelTransform(N=50)
    res = ht.transform(lambda x: 1. / x, k, False)
    assert len(res) == 10


def test_k_scalar():
    k = 1
    ht = HankelTransform(N=50)
    res = ht.transform(lambda x: 1. / x, k, False)
    assert np.isscalar(res)


def test_ret_err():
    ht = HankelTransform(N=50)
    res, err = ht.integrate(lambda x: 1. / x, True)


def test_ret_cumsum():
    ht = HankelTransform(N=50)
    res, cumsum = ht.integrate(lambda x: 1. / x, False, True)


def test_ret_err_and_cumsum():
    ht = HankelTransform(N=50)
    res, err, cumsum = ht.integrate(lambda x: 1. / x, True, True)


def test_equivalence_of_integrate_and_transform():
    ht = HankelTransform(N=50)
    int = ht.integrate(lambda x: 1, False)
    tr = ht.transform(lambda x: 1. / x, ret_err=False)
    assert int == tr
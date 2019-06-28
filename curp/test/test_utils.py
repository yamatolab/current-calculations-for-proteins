import numpy
from nose.tools import *
import numpy.testing as npt

from curp import utility


def gen_tensor3(xs, ys):
    for x, y in zip(xs, ys):
        yield numpy.array([[x[0]*y[0], x[0]*y[1], x[0]*y[2]],
                           [x[1]*y[0], x[1]*y[1], x[1]*y[2]],
                           [x[2]*y[0], x[2]*y[1], x[2]*y[2]]])

def test_tensor1_same():
    xs = numpy.array([[0.11,0.08,-0.02]])
    res = utility.tensor(xs, xs)
    right = list(gen_tensor3(xs, xs))

    npt.assert_almost_equal(res, right)

def test_tensor1_differ():
    xs = numpy.array([[0.11,0.08,-0.02]])
    ys = numpy.array([[0.21,-0.10,0.05]])
    res = utility.tensor(xs, ys)
    right = list(gen_tensor3(xs, ys))

    npt.assert_almost_equal(res, right)


def test_tensor3_same():
    xs = numpy.array([[0.11,0.08,-0.02],
                      [0.21,-0.10,0.05],
                      [0.05,-0.01,-0.13]])
    res = utility.tensor(xs, xs)
    right = list(gen_tensor3(xs, xs))

    npt.assert_almost_equal(res, right)

def test_tensor3_differ():
    xs = numpy.array([[0.11,0.08,-0.02],
                      [0.21,-0.10,0.05],
                      [0.05,-0.01,-0.13]])
    ys = numpy.array([[0.09,0.15,-0.22],
                      [-0.21,-0.19,0.15],
                      [-0.04,-0.08,-0.12]])
    res = utility.tensor(xs, ys)
    right = list(gen_tensor3(xs, ys))

    npt.assert_almost_equal(res, right)


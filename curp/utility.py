
# import numpy
# class DimensionError(Exception): pass
# class RangeError(Exception): pass
# class MultiArray(dict):

#     def __init__(self, dim=2, default=0.0):
#         self.__ndim     = dim
#         self.__default = default

#     def __getitem__(self, key):
#         self._check_key(key)
#         return self.get(key, self.__default)

#     def __setitem__(self, key, value):
#         self._check_key(key)
#         super(self.__class__, self).__setitem__(key, value)

#     def _check_key(self, key):
#         if len(key) != self.__ndim:
#             raise DimensionError(key)
#         for k in key:
#             if not isinstance(k, int):
#                 raise TypeError(k)
#             if k < 0:
#                 raise RangeError(k)

#     def array(self):
#         dims = zip(*self.keys())
#         dim_maxs = [max(dim) for dim in dims]
#         # dim_mins = [min(dim) for dim in dims]

#         array = numpy.zeros(dim_maxs, dtype=float)
#         for dmax in dim_maxs:

import numpy

def tensor(xs, ys):
    axes_x = xs.shape
    axes_y = ys.shape 

    if axes_x != axes_y:
        print(axes_x, axes_y)
        raise ValueError

    dim0 = axes_x[-1]

    if len(axes_x) == 1:
        result = numpy.zeros([dim0, dim0])
    elif len(axes_x) == 2:
        result = numpy.zeros([axes_x[0], dim0, dim0])
    elif len(axes_x) == 3:
        result = numpy.zeros([axes_x[0], axes_x[1], dim0, dim0])
    elif len(axes_x) == 4:
        result = numpy.zeros([axes_x[0], axes_x[1], axes_x[2], dim0, dim0])
    else:
        raise numpy.DimensionError()

    for i in range(dim0):
        for j in range(dim0):
            result[..., i, j] = xs[..., i] * ys[..., j]

    return result

def tensor1(x, y):
    result = numpy.zeros([3,3])
    for i in range(3):
        for j in range(3):
            result[i, j] = x[i] * y[j]

    return result

def tensor2(xs, ys):
    axes_x = xs.shape
    axes_y = ys.shape 

    if axes_x != axes_y:
        raise ValueError

    dim0 = axes_x[-1]
    if len(axes_x) == 2:
        result = numpy.zeros([axes_x[0], dim0, dim0])
        for n, (x, y) in enumerate(zip(xs,ys)):
            for i in range(dim0):
                for j in range(dim0):
                    result[n, i, j] = x[i] * y[j]

        return result


def time_tensor(natom=10**5):
    xs = numpy.ones([natom,3])
    ys = numpy.ones([natom,3])

    from benchmarker import Benchmarker
    with Benchmarker(width=20) as bm:

        with bm('tensor'):
            result = tensor2(xs, ys)

        with bm('tensor 1'):
            result = numpy.zeros([natom,3,3])
            for i,(x, y) in enumerate(zip(xs, ys)):
                result[i] = tensor1(x,y)
            # print(result)

        with bm('tensor2'):
            result = tensor3(xs, ys)

def time_tensor2(natom=10**5):
    xs = numpy.ones([100, natom,3])
    ys = numpy.ones([100, natom,3])
    # xs = numpy.array([1,2,3])
    # ys = numpy.array([4,5,6])

    from benchmarker import Benchmarker
    with Benchmarker(width=10) as bm:
        with bm('tensor'):
            result = tensor(xs, ys)

        with bm('dot product'):
            aa = numpy.arange(10**5)
            bb = numpy.arange(10**5)
            aabb = numpy.dot(aa, bb)
        print(result)
        print(aabb)

# def tensor(xs, ys):
#     result = numpy.zeros(xs.shape)
#     for x, y in zip(xs, ys):
#          = x[0]*y[0], x[0]*y[1], x[0]*y[2]

# from scipy import sparse
class InvalidShapeError: pass
class TBFMatrix_Old:
    """
    The old class for two-body force
    """

    def __init__(self, shape):
        self.__shape = shape
        self.xs = sparse.lil_matrix(shape)
        self.ys = sparse.lil_matrix(shape)
        self.zs = sparse.lil_matrix(shape)

    def get_shape(self):
        return self.__shape

    def __setitem__(self, index, value):
        self.xs[index] = value[0]
        self.ys[index] = value[1]
        self.zs[index] = value[2]

    def __getitem__(self, index):
        return numpy.array([self.xs[index],self.ys[index],self.zs[index]])

    def __add__(self, other):
        if self.__shape != other.get_shape():
            raise InvalidShapeError
        new = TBFMatrix(self.__shape)
        new.xs = self.xs + other.xs
        new.ys = self.ys + other.ys
        new.zs = self.zs + other.zs
        return new

    def iadd(self, index, value):
        self.xs[index] += value[0]
        self.ys[index] += value[1]
        self.zs[index] += value[2]

    def isub(self, index, value):
        self.xs[index] -= value[0]
        self.ys[index] -= value[1]
        self.zs[index] -= value[2]

    def __iadd__(self, other):
        if self.__shape != other.get_shape():
            raise InvalidShapeError
        self.xs += other.xs
        self.ys += other.ys
        self.zs += other.zs
        return self

    def minus(self, value):
        self.xs[index] -= value[0]
        self.ys[index] -= value[1]
        self.zs[index] -= value[2]

    def __isub__(self, value):
        self.xs[index] -= value[0]
        self.ys[index] -= value[1]
        self.zs[index] -= value[2]

    def items(self):
        ilist, jlist = self.xs.nonzero()
        for i, j in zip(ilist, jlist):
            yield i, j, numpy.array([
                self.xs[i,j], self.ys[i,j], self.zs[i,j] ])

    def __repr__(self):
        mes = ( '[{i},{j}] = {value}'.format(i=i,j=j,value=str(ary))
                for i, j, ary in self.items() )

        return '\n'.join(mes)


class TBFMatrix(dict):

    def __init__(self, shape):
        dict.__init__(self)
        self.__shape = shape

    def get_shape(self):
        return self.__shape

    def __getitem__(self, key):
        return self.get(key, numpy.zeros([3]))

    def __add__(self, other):
        if self.__shape != other.get_shape():
            raise InvalidShapeError
        new = self.copy()
        for key, value in other.items():
            new[key] += value
        return new

    def tadd(self, key, value):
        i, j = key
        self[i, j] += value
        self[j, i] -= value


import functools
def timeit(fun):
    import time
    @functools.wraps(fun)
    def wrapped(*args, **kwargs):
        t0 = time.time()
        result = fun(*args, **kwargs)
        t1 = time.time()
        print("{0} : {1}").format(fn.__name__, t1 - t0)
        return result
    return wrapped

class TimeStore:

    def __init__(self):
        self.__tag_to_times = {}
        self.__tags = []

    def store_time(self, tag, time):
        if tag not in self.__tags:
            self.__tag_to_times[tag] = []
            self.__tags.append(tag)

        self.__tag_to_times[tag].append(time)

    def gen_time_info(self):
        for tag in self.__tags:
            yield tag, self.__tag_to_times[tag]


if __name__ == '__main__':
    time_tensor2()

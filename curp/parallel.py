"""parallel - Classes to handle processor, parallel or serial

Classes:
ParallelProcessor -- Run the calculations in parallel.
SequentialProcess -- Run the calculations in serial.
"""

from __future__ import print_function

import itertools as it
import time

try:
    from mpi4py import MPI
    use_mpi = True
except ImportError:
    use_mpi = False

################################################################################
class ParallelProcessor:

    """
    The class to run in parallel from given method and dictionary list,
    using mpi4py library.

    >>> para = ParallelProcessor()
    >>> results = para.run(method, crd=crds, vel=vels)

    The calculated result return as generator.

    >>> for res in results:
        if para.is_root():
    ...     print(res)

    To perform test calculation, You can be run as the main module:

        $ mpiexec -n 8 python parallel.py

    """

    def __init__(self, chunksize=5):
        self.__chunksize = chunksize

        self.__comm = MPI.COMM_WORLD
        self.__nproc = self.__comm.Get_size()
        self.__iproc = self.__comm.Get_rank()

        self.__nstep = 0
        self.__other_time    = None
        self.__scatter_times = []
        self.__body_times    = []
        self.__gather_times  = []

    def is_root(self):
        return self.__iproc == 0

    def get_rank(self):
        return self.__iproc

    def get_nproc(self):
        return self.__nproc

    def _gen_kwds(self, kwds):
        """Generate each dictionary from listed dictionaries."""

        ks = kwds.keys()    # kwd1, kwd2, ...
        npair = len(ks)     # the number of kwds

        if self.is_root():
            vs = kwds.values()  # val1_iter, val2_iter, ...
            pair_vs_iter = it.izip(*vs)

        else:
            def gen_none():
                while True:
                    yield npair*[None]

            pair_vs_iter = gen_none()

        while True:

            is_end = False
            if self.is_root():
                try:
                    pair_vs = pair_vs_iter.next()
                except StopIteration:
                    is_end = True

            else:
                pair_vs = pair_vs_iter.next()

            is_end = self.__comm.bcast(is_end, root=0)
            
            if is_end:
                break

            # (val1, val2, ...)_1, (val1, val2, ...)_2, ...
            yield dict( (kwd, val) for kwd, val in zip(ks, pair_vs) )

    def _gen_kwds_nproc_list(self, kwds_iter):
        """Generate a list of the number of precesses for keyword dictionary."""
        iproc = 0
        kwds_list_sumup = []
        for kwds in kwds_iter:
            kwds_list_sumup.append(kwds)
            iproc += 1

            if iproc == self.__nproc:
                yield kwds_list_sumup
                kwds_list_sumup = []
                iproc = 0

        else:
            if len(kwds_list_sumup) != 0:
                for ip in range(self.__nproc-len(kwds_list_sumup)):
                    kwds_list_sumup.append(None)
                yield kwds_list_sumup

    def run(self, method, *args, **kwds):
        """Generate result after run in parallell."""

        scatter_times = []
        method_times  = []
        gather_times  = []
        other_time = 0.0

        kwds_nproc_iter = self._gen_kwds_nproc_list( self._gen_kwds(kwds) )

        t0 = time.time()

        istep = 0
        for istep, kwds_list in enumerate(kwds_nproc_iter):

            # scatter method's kwds on each process

            t1 = time.time()
            kwds = self.__comm.scatter(kwds_list, root=0)

            # run on each process
            t2 = time.time()
            if kwds:
                result = method(**kwds)
            else:
                result = None

            # gather return data from on each process
            t3 = time.time()
            if self.get_nproc() != 1:
                results = self.__comm.gather(result, root=0)
            else:
                results = [result]
            t4 = time.time()

            # split all results on each processes.
            if results:
                for res in results:
                    if res is not None:
                        yield res

            # store time intervals for each process
            scatter_times.append(t2-t1)
            method_times.append(t3-t2)
            gather_times.append(t4-t3)

            t5 = time.time()
            other_time += t5 - t4

        t6 = time.time()

        self.__mpi_time   = t6 - t0
        self.__other_time = other_time
        self.__nstep = istep
        self.__scatter_times = self.__comm.gather(scatter_times, root=0)
        self.__body_times  = self.__comm.gather(method_times, root=0)
        self.__gather_times  = self.__comm.gather(gather_times, root=0)

    def run_nogather(self, method, *args, **kwds):
        """Generate result after run in parallell."""

        scatter_times = []
        method_times  = []
        gather_times  = []

        kwds_nproc_iter = self._gen_kwds_nproc_list( self._gen_kwds(kwds) )

        t0 = time.time()

        istep = 0
        for istep, kwds_list in enumerate(kwds_nproc_iter):

            # scatter method's kwds on each process

            t1 = time.time()
            kwds = self.__comm.scatter(kwds_list, root=0)

            # run on each process
            t2 = time.time()
            if kwds:
                method(**kwds)

            # store time intervals for each process
            scatter_times.append(t2-t1)
            method_times.append(t3-t2)
            gather_times.append(0.0)

            # yield result
            yield None

        t5 = time.time()
        self.__mpi_time = t5 - t0
        self.__other_time = 0.0
        self.__nstep = istep
        self.__scatter_times = self.__comm.gather(scatter_times, root=0)
        self.__body_times  = self.__comm.gather(method_times, root=0)
        self.__gather_times  = self.__comm.gather(gather_times, root=0)

    def get_time_info(self, nstep):

        import numpy
        nproc = self.__nproc
        denom = 1.0/float(nproc*nstep)

        stime_total = sum([sum(times) for times in self.__scatter_times])
        btime_total = sum([sum(times) for times in self.__body_times   ])
        gtime_total = sum([sum(times) for times in self.__gather_times ])

        return [ ('Scatter', stime_total*denom),
                 ('Body'   , btime_total*denom),
                 ('Gather' , gtime_total*denom) ]

    def get_mpi_time(self):
        return self.__mpi_time

    def get_other_time(self):
        return self.__other_time

    def write(self, *args, **kwds):
        if self.is_root():
            print(*args, **kwds)

################################################################################
class SequentialProcessor:

    """Class to run in serial from given method and dictionary list,
    when there is not mpi4py library.

    >>> para = SequentialProcessor()
    >>> results = para.run(method, crd=crds, vel=vels)

    The calculated result return as generator.

    >>> for res in results:
        if para.is_root():
    ...     print(res)

    To perform test calculation, You can be run as the main module:

        $ python parallel.py

    """

    def __init__(self):
        self.__iproc = 0
        self.__body_times = []

    def is_root(self):
        return True

    def get_rank(self):
        return 0

    def get_nproc(self):
        return 1

    def _gen_kwds(self, kwds):
        """Generate each dictionary from listed dictionaries."""

        ks = kwds.keys()    # kwd1, kwd2, ...
        npair = len(ks)     # the number of kwds
        vs = kwds.values()  # val1_iter, val2_iter, ...
        pair_vs_iter = it.izip(*vs)

        for pair_vs in pair_vs_iter:
            # (val1, val2, ...)_1, (val1, val2, ...)_2, ...
            yield dict( (ks[i], pair_vs[i]) for i in range(npair) )

    def run(self, method, *args, **kwds):
        """Generate result after run in Serial."""

        method_times = []

        kwds_iter = self._gen_kwds(kwds)
        istep = None
        for istep, key_to_datas in enumerate(self._gen_kwds(kwds)):

            t1 = time.time()
            result = method(**key_to_datas)
            t2 = time.time()

            # store time intervals for each process
            method_times.append(t2-t1)

            yield result

        self.__nstep = istep if istep else 0
        self.__body_times  = method_times

    def get_time_info(self, nstep):
        t = sum(time for time in self.__body_times) / float(nstep)
        return [ ('Body' , t ) ]

    def get_other_time(self):
        return None

    def write(self, *args, **kwds):
        print(*args, **kwds)


def main():


    def method(crd, vel):
        import time
        time.sleep(0.1)
        return crd*vel

    import time
    t1 = time.time()

    crds = ( i+1 for i in range(1000))
    vels = ( 10*(i+1) for i in range(1000))
    para = ParallelProcessor()

    if para.is_root():
        t2 = time.time()
        print('elasped time : t2 - t1 : {}'.format(t2-t1))

    results = para.run(method, crd=crds, vel=vels)
    if para.is_root():
        t3 = time.time()
        print('elasped time : t3 - t2 : {}'.format(t3-t2))

    for r in results:
        # if para.is_root():
        r

    if para.is_root():
        t4 = time.time()
        print('elasped time : t4 - t3 : {}'.format(t4-t3))

def main_list(method, natom, ntraj):

    import time
    t1 = time.time()

    # initialization
    import copy
    crd  = [ (iatm+1, iatm+2, iatm+3) for iatm in range(natom) ]
    crds = ( copy.copy(crd) for itraj in range(ntraj) )
    vel  = [ (10*(iatm+1), 10*(iatm+2), 10*(iatm+3)) for iatm in range(natom) ]
    vels = ( copy.copy(crd) for itraj in range(ntraj) )

    run_bench(method, crds, vels, 'List', t1)

def main_array(method, natom, ntraj):

    import time
    t1 = time.time()

    # initialization
    import copy
    import numpy
    crd = numpy.ones([natom,3])

    crds = ( crd.copy() for itraj in range(ntraj) )
    vel = numpy.ones([natom,3])
    vels = ( vel.copy() for itraj in range(ntraj) )

    run_bench(method, crds, vels, 'Array', t1)

def main_arraygen(natom, ntraj):
    # initialization
    import copy
    import numpy

    par = ParallelProcessor()
    label = 'Array 2'

    import time
    t1 = time.time()

    def method(data):
        import time
        time.sleep(0.5)
        itraj, (crd, vel) = data
        return itraj

    if par.is_root():

        crdvel = numpy.ones([natom,3])

        def gen_data_iter(ntraj):
            for itraj in range(ntraj):
                yield itraj, (crdvel.copy(), crdvel.copy())

        data_iter = gen_data_iter(ntraj)
    else:
        data_iter = None
        # def gen_none_iter():
            # while True:
                # yield None
        # data_iter = gen_none_iter()

        # def gen_none_iter(ntraj):
            # for itraj in range(ntraj):
                # yield None

        # data_iter = gen_none_iter(ntraj)

    t2 = time.time()
    msg = 'elasped time : {} preparation     : {:8.5f}'.format(label, t2-t1)
    par.write(msg)

    gen_result = par.run(method, data=data_iter)

    t3 = time.time()
    msg = 'elasped time : {} passed into run : {:8.5f}'.format(label, t3-t2)
    par.write(msg)

    # if par.is_root():
    for r in gen_result:
        # par.write(r) 
        # if par.is_root():
            print(r)

    t4 = time.time()
    msg = 'elasped time : {} perform the run : {:8.5f}'.format(label, t4-t3)
    par.write(msg)


def main_arraygen(natom, ntraj):
    # initialization
    import copy
    import numpy

    par = ParallelProcessor()
    label = 'Array 2'

    import time
    t1 = time.time()

    def method(data):
        import time
        time.sleep(0.5)
        itraj, (crd, vel) = data
        return itraj

    if par.is_root():

        crdvel = numpy.ones([natom,3])

        def gen_data_iter(ntraj):
            for itraj in range(ntraj):
                yield itraj, (crdvel.copy(), crdvel.copy())

        data_iter = gen_data_iter(ntraj)
    else:
        data_iter = None
        # def gen_none_iter():
            # while True:
                # yield None
        # data_iter = gen_none_iter()

        # def gen_none_iter(ntraj):
            # for itraj in range(ntraj):
                # yield None

        # data_iter = gen_none_iter(ntraj)

    t2 = time.time()
    msg = 'elasped time : {} preparation     : {:8.5f}'.format(label, t2-t1)
    par.write(msg)

    gen_result = par.run(method, data=data_iter)

    t3 = time.time()
    msg = 'elasped time : {} passed into run : {:8.5f}'.format(label, t3-t2)
    par.write(msg)

    # if par.is_root():
    for r in gen_result:
        # par.write(r) 
        # if par.is_root():
            print(r)

    t4 = time.time()
    msg = 'elasped time : {} perform the run : {:8.5f}'.format(label, t4-t3)
    par.write(msg)


def run_bench(method, crds, vels, label, t):

    import time

    para = ParallelProcessor()

    t2 = time.time()
    msg = 'elasped time : {} preparation     : {:8.5f}'.format(label, t2-t)
    para.write(msg)

    gen_result = para.run(method, crd=crds, vel=vels)
    t3 = time.time()
    msg = 'elasped time : {} passed into run : {:8.5f}'.format(label, t3-t2)
    para.write(msg)

    for r in gen_result:
        # if para.is_root():
        # para.write(r)
        para.write(para.get_rank(),r)

    t4 = time.time()
    msg = 'elasped time : {} perform the run : {:8.5f}'.format(label, t4-t3)
    para.write(msg)


if __name__ == '__main__':
    natom = 10**5
    ntraj = 16

    def method(crd, vel):
        import time
        time.sleep(0.5)
        # import utility
        # return utility.tensor(crd, vel)
        return len(crd)


    main_list(method, natom, ntraj)
    main_array(method, natom, ntraj)
    main_arraygen(natom, ntraj)

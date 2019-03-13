from __future__ import print_function

import itertools

################################################################################
class InteractionTableGenerator:

    """

    data structure:
        [jint, (iatm, jatm_begin, jatm_end)]

    >>> table1 = InteractionTable(natom)
    >>> table1
    <InteractionTable>

    exclude the bonded pair
    >>> table2 = table1.filter(if_excluded_bond)
    <InteractionTable>

    >>> list(table2)

    perform cutoff with cutoff length to make a interaction table within
    length from each atoms
    >>> table3 = table2.cutoff(crd, 10.0)
    >>> table.map()
    """

    def __init__(self, natom=None, base_table=None):
        if natom:
            self.__base_generater = self._init_table(natom)
        else:
            if base_table:
                self.__base_generater = base_table
            else:
                self.__base_generater = []

        self.__filters = []
        self.__cache = None

    def _init_table(self, natom):
        for iatm_1 in range(natom-1):
            iatm = iatm_1 + 1
            yield (iatm, iatm+1, natom)

    def filter(self, if_fun):
        self.__filters.append(if_fun)
        return self 

    def insert_filter(self, if_fun):
        self.__filters.insert(0, if_fun)
        return self

    def list(self):
        return list(self.__apply())

    def __iter__(self):
        return self

    def next(self):
        return self.__apply().next()
    __next__ = next # version 3.x

    def copy(self):
        gen1, gen2 = itertools.tee(self.__apply())
        self.__apply = gen1
        self.__filters = []
        return InteractionTable(base_table=gen2)

    def __apply(self):
        """Generator: apply if_function to the table."""

        if self.__filters == []:
            for iatm, jatm_beg, jatm_end in self.__base_generater:
                yield iatm, jatm_beg, jatm_end

        for iatm, jatm_beg, jatm_end in self.__base_generater:

            jatm_beg_new = 0
            jatm_end_new = jatm_beg
            
            for jatm in range(jatm_beg, jatm_end+1):
                print('jatm before', jatm, jatm_end)

                for if_fun in self.__filters:
                    if not if_fun(iatm, jatm):
                        print('aaa')
                        if jatm_beg_new > 0:
                            print('bbb',iatm, jatm_beg_new, jatm_end_new)
                ########################### here bug ########################
                            yield iatm, jatm_beg_new, jatm_end_new
                            print('ccc')
                            jatm_beg_new = 0
                        print('ddd')
                        # break

                else: # filters are all ok.
                    print('filter else', jatm, jatm_end)
                    if jatm_beg_new==0: jatm_beg_new = jatm
                    jatm_end_new = jatm

                print('jatm after', jatm)

            else:
                print('jatm else', jatm)
                if jatm_beg_new > 0:
                    yield iatm, jatm_beg_new, jatm_end_new

    def cutoff(self, crd, length):
        return self.filter(self._if_cutoff(crd, length))

    def __filter_directly(self, if_function):
        """Apply if_function to the table, generator version."""

        for iatm, jatm_beg, jatm_end in self.__base_generater:

            jatm_beg_new = 0
            jatm_end_new = jatm_beg
            
            for jatm in range(jatm_beg, jatm_end+1):

                if if_function(iatm, jatm):
                    if jatm_beg_new==0: jatm_beg_new = jatm
                    jatm_end_new = jatm
                else:
                    if jatm_beg_new==0: continue
                    yield iatm, jatm_beg_new, jatm_end_new
                    jatm_beg_new = 0

            else:
                if jatm_beg_new > 0:
                    yield iatm, jatm_beg_new, jatm_end_new

    def _if_cutoff(self, crd, length):
        length2 = length**2

        def wrapper(iatm, jatm):
            dist2 = numpy.dot(crd[iatm,:] - crd[jatm,:])
            return dist2 < length2

        return wrapper

    def gen_each_atom(self):
        """Generator: InteractionTable data divided by each atoms."""
        cur_iatm  = 1
        int_table = []

        for iatm, jatm_beg, jatm_end in self:

            if iatm == cur_iatm:
                int_table.append((iatm, jatm_beg, jatm_end))
            else:
                yield int_table
                cur_iatm += 1

                while cur_iatm < iatm:
                    yield []
                    cur_iatm += 1

                int_table = [(iatm, jatm_beg, jatm_end)]

        else:
            yield int_table

    def get_each_atom(self):
        if self.__cache is None:
            return list(self.gen_each_atom())
        else:
            return self.__cache


################################################################################
class InteractionTableList:

    """

    data structure:
        [jint, (iatm, jatm_begin, jatm_end)]

    >>> table1 = InteractionTable(natom)
    >>> table1
    <InteractionTable>

    exclude the bonded pair
    >>> table2 = table1.filter(if_excluded_bond)
    <InteractionTable>

    >>> list(table2)

    perform cutoff with cutoff length to make a interaction table within
    length from each atoms
    >>> table3 = table2.cutoff(crd, 10.0)
    >>> table.map()
    """

    _memory_limit = 10 # MB

    def __init__(self, natom=None, base_table=None):
        self.__natom = natom
        if self.__natom:
            self.__base_generater = self._init_table(self.__natom)
        else:
            if base_table:
                self.__base_generater = base_table
            else:
                self.__base_generater = []

        self.__filters = []
        self.__ifilters = []
        self.__cnt = -1
        self.__num_interacts = 0
        self.__table = []

    def _init_table(self, natom):
        for iatm_1 in range(natom-1):
            iatm = iatm_1 + 1
            yield (iatm, iatm+1, natom)

    def filter_iatom(self, if_fun):
        self.__ifilters.append(if_fun)
        return self

    def filter(self, if_fun):
        self.__filters.append(if_fun)
        return self 

    def filter_or(self, iatm_beg, iatm_end):
        self.__filters.append( ('or', iatm_beg, iatm_end) )
        return self

    def filter_and(self, iatm_beg, iatm_end):
        self.__filters.append( ('and', iatm_beg, iatm_end) )
        return self

    def filter_fun(self, if_fun):
        self.__filters.append( ('fun', if_fun) )
        return self

    def insert_filter(self, if_fun):
        self.__filters.insert(0, if_fun)
        return self

    def list(self):
        return self.__apply()

    def __iter__(self):
        return self

    def next(self):
        self.__cnt += 1

        if self.__num_interacts == 0:
            self.__apply()

        if self.__cnt >= self.__num_interacts:
            raise StopIteration()

        return self.__table[self.__cnt]

    __next__ = next # version 3.x

    def copy(self):
        gen1, gen2 = itertools.tee(self.__apply())
        self.__apply = gen1
        self.__filters = []
        self.__ifilters = []
        return InteractionTable(base_table=gen2)

    def __apply(self):
        """Apply if_function to the table."""

        def gen_filter_or(table, i_beg, i_end):
            for iatm, jatm_beg, jatm_end in table:

                if iatm < i_beg:
                    yield iatm, i_beg, i_end

                elif i_beg <= iatm <= i_end:
                    yield iatm, jatm_beg, jatm_end

                elif i_end < iatm:
                    break 

                else:
                    pass

        def gen_filter_and(table, i_beg, i_end):
            for iatm, jatm_beg, jatm_end in table:

                if iatm < i_beg:
                    pass

                elif i_beg <= iatm <= i_end-1:
                    yield iatm, jatm_beg, jatm_end

                elif i_end-1 < iatm:
                    break 

                else:
                    pass

        def gen_filter_fun(table, if_fun):

            for iatm, jatm_beg, jatm_end in table:
                
                jatm_beg_new = 0
                jatm_end_new = jatm_beg
                
                for jatm in range(jatm_beg, jatm_end+1):

                    if not if_fun(iatm, jatm):
                        if jatm_beg_new > 0:
                            yield iatm, jatm_beg_new, jatm_end_new
                            jatm_beg_new = 0

                    else:
                        if jatm_beg_new==0: jatm_beg_new = jatm
                        jatm_end_new = jatm

                else:
                    if jatm_beg_new > 0:
                        yield iatm, jatm_beg_new, jatm_end_new

        table = self.__base_generater
        for f in self.__filters:
            ftype = f[0]

            if ftype == 'or':
                i_beg, i_end = f[1], f[2]
                next_table = gen_filter_or(table, i_beg, i_end)

            elif ftype == 'and':
                i_beg, i_end = f[1], f[2]
                next_table = gen_filter_and(table, i_beg, i_end)

            elif ftype == 'fun':
                if_fun = f[1]
                next_table = gen_filter_fun(table, if_fun)

            else:
                pass

            table = next_table

        self.__table = list(table)
        self.__num_interacts = len(self.__table)

    def cutoff(self, crd, length):
        return self.filter(self._if_cutoff(crd, length))

    def _if_cutoff(self, crd, length):
        length2 = length**2

        def wrapper(iatm, jatm):
            dist2 = numpy.dot(crd[iatm,:] - crd[jatm,:])
            return dist2 < length2

        return wrapper

    def get_size_mem(self, iatm_min, iatm_max, jatm_min, jatm_max):
        scale = 1
        ibound = iatm_max - iatm_min + 1
        jbound = jatm_max - jatm_min + 1
        memory = 8 * ibound * jbound * 3 * scale
        return memory

    def get_size_pair(self, npair):
        scale = 1
        memory = 8 * npair * 3 * scale
        return memory

    get_size = get_size_pair

    def gen_each_atom(self):
        """Generator: InteractionTable data divided by each atoms."""
        cur_iatm  = 1
        table = []

        for iatm, jatm_beg, jatm_end in self:

            if iatm == cur_iatm:
                table.append((iatm, jatm_beg, jatm_end))
            else:
                yield table
                cur_iatm += 1

                while cur_iatm < iatm:
                    yield []
                    cur_iatm += 1

                table = [(iatm, jatm_beg, jatm_end)]

        else:
            yield table

    def gen_divided_table_old(self, num_pair_limits=None):
        """Generate a table divided to each number of pair."""

        if num_pair_limits is None:
            num_pair_limits = self._limit_pair

        npair = 0
        table = []

        for iatm, jatm_beg, jatm_end in self:
            npair += jatm_end - jatm_beg + 1
            table.append( (iatm, jatm_beg, jatm_end) )

            if npair > num_pair_limits:
                yield table
                npair = 0
                table = []

        else:
            if table:
                yield table

    def gen_divided_table_pair(self, mem_limit=None):
        """Generate a table divided to each number of pair."""

        self.__memories = []

        if mem_limit is None:
            mem_limit = self._memory_limit

        mem_lower = mem_limit * 10**6
        mem_upper = 1.2 * mem_lower

        npair = 0
        table = []
        npair_p = 0

        base_table = list(self)

        for itab, (iatm, jatm_beg, jatm_end) in enumerate(base_table):

            table += [(iatm, jatm_beg, jatm_end)]
            npair += jatm_end - jatm_beg + 1
            mem = self.get_size(npair)

            if mem_lower <= mem <= mem_upper:
                self.__memories.append( mem )
                yield table
                npair = 0
                table = []

            elif mem_upper < mem:
                iatm_p, jatm_beg_p, jatm_end_p = table.pop()
                mem = self.get_size(npair_p)
                self.__memories.append( mem )
                yield table
                npair = jatm_beg_p - jatm_beg_p + 1
                table = [(iatm_p, jatm_beg_p, jatm_end_p)]

            npair_p = npair

        else:
            mem = self.get_size(npair_p)
            self.__memories.append( mem )
            if table: yield table

    def gen_divided_table_memory(self):
        """Generate a table divided by memory limit."""

        # if self._memory_limit is None:
            # memory_limt = self._memory_limit

        self.__memories = []
        memory_lower = self._memory_limit * 10**3
        memory_upper = 1.2 * memory_lower
        base_table = list(self)

        for itab, (iatm,jatm_beg,jatm_end) in enumerate(base_table):

            if itab==0:
                iatm_min, iatm_max = iatm, iatm
                jatm_min, jatm_max = jatm_beg, jatm_end
                table = [(iatm, jatm_beg, jatm_end)]
                memory = 0
                continue

            table += [(iatm, jatm_beg, jatm_end)]
            iatm_max = iatm
            jatm_min = min(jatm_beg, jatm_min)
            jatm_max = max(jatm_end, jatm_max)

            memory = self.get_size(iatm_min, iatm_max, jatm_min, jatm_max)

            if memory_lower < memory < memory_upper:
                self.__memories.append( memory )
                # print('memory', iatm_min, iatm_max, jatm_min, jatm_max)
                print('bound', iatm_max-iatm_min+1, jatm_max-jatm_min+1)
                yield table
                iatm_min = iatm
                jatm_min = jatm_end
                table = []

            elif memory_upper <= memory:
                # print('memory', iatm_min, iatm_max, jatm_min, jatm_max)
                iatm_p, jatm_beg_p, jatm_end_p  = table.pop()
                # print('bound', iatm_max-iatm_min+1, jatm_max-jatm_min+1)
                memory = self.get_size(
                        iatm_min_p, iatm_max_p, jatm_min_p, jatm_max_p)
                print('bound', iatm_max_p-iatm_min_p+1, jatm_max_p-jatm_min_p+1)
                self.__memories.append( memory )
                yield table
                iatm_min = iatm_p
                jatm_min = jatm_end_p
                table = [(iatm_p, jatm_beg_p, jatm_end_p)]

            iatm_min_p = iatm_min
            iatm_max_p = iatm_max
            jatm_min_p = jatm_min
            jatm_max_p = jatm_max
            
        else:
            memory = self.get_size(iatm_min, iatm_max, jatm_min, jatm_max)
            self.__memories.append( memory )
            if table: yield table

    gen_divided_table = gen_divided_table_pair

    def get_table_memories(self):
        return self.__memories


InteractionTable = InteractionTableList

################################################################################
def make_accesible_table(pairs, natom):
    """Convert the given pair list to easily accessible table.
    pair_table: iatm => (jatm1, jatm2, ...)
    """
    # recreate the bonded pair list to make non-bonded interaction table
    # effectively.
    pair_table = [ [] for iatm in range(natom) ]
    for iatm, jatm in pairs:
        pair_table[iatm-1].append(jatm)
    return pair_table


if __name__ == '__main__':

    target = [5,7,8,9]

    def if_fun_or(iatm, jatm):
        return (iatm in target) or (jatm in target)

    def if_fun_and(iatm, jatm):
        return (iatm in target) and (jatm in target)

    # def if_fun2(iatm, jatm):
    #     return 101 <= iatm <= 200

    # def if_fun3(iatm):
    #     return iatm in range(100, 200)


    natom = 5000
    from benchmarker import Benchmarker
    with Benchmarker(width=20) as bm:

        with bm('or'):
            print()
            table1 = InteractionTable(natom)
            t = table1.filter_or(5, 10)
            for a in t.gen_divided_table(5000):
                print(a)

        with bm('and'):
            print()
            table2 = InteractionTable(natom)
            t = table2.filter_and(10, 20)
            for a in t.gen_divided_table(40):
                print(a)
                
        # with bm('if_fun_or'):
        #     print()
        #     table3 = InteractionTable(natom)
        #     t = table3.filter_fun(if_fun_or)
        #     for a in t.gen_divided_table(40):
        #         print(a)

        # with bm('if_fun_and'):
        #     print()
        #     table4 = InteractionTable(natom)
        #     t = table4.filter_fun(if_fun_and)
        #     for a in t.gen_divided_table(40):
        #         print(a)

        with bm('or, if_fun_or'):
            print()
            table = InteractionTable(natom)
            t = table.filter_or(5, 10).filter_fun(if_fun_or)
            for a in t.gen_divided_table(40):
                print(a)


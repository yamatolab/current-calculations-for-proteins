#! /usr/bin/env python
from __future__ import print_function
import sys
import math
import numpy


class TensorParser:

    section_end = '%end'

    def __init__(self, filename):
        self.__filename = filename

    def parse(self):
        with open(self.__filename, 'rb') as file:
            names, tensors = self.parse_content(file)
        return names, tensors

    def parse_content(self, file):

        tensors = []

        gen_lines = self.gen_optimized_lines(file)
        sections = []
        for line in gen_lines:
            if not self.is_section(line): continue

            flagname = line[1:].strip()
            if flagname == 'data':
                lines = list( self.gen_section_lines(gen_lines) )
                names, tensor = self.parse_data(lines)
                tensors.append( tensor )

        return names, tensors

    def parse_data(self, lines):
        names = []
        tensors = []
        for line in lines:
            cols = line.split()
            names.append(cols[0])
            v = [float(c) for c in cols[1:]]
            tensor = [ [v[0], v[1], v[2]],
                       [v[3], v[4], v[5]],
                       [v[6], v[7], v[8]] ]
            tensors.append(tensor)

        # array = numpy.array(v)
        # print(array.reshape(3,3,len(array)))
        return names, numpy.array(tensors)

    def gen_optimized_lines(self, file):
        # first flag
        for line in file:
            if self.is_section(line):
                yield line
                break
            yield line

        for line in file:
            if self.is_section(line):
                yield self.section_end
                yield line
            yield line

        else:
            yield self.section_end

    def split_content(self, gen):
        lines = []
        for line in gen:
            if line.startswith('%'):
                pass

    def is_section(self, line):
        if line == self.section_end:
            return False
        elif line.startswith('%'):
            return True
        else:
            return False

    def gen_section_lines(self, gen):
        for line in gen:
            if line == self.section_end:
                break
            elif self.is_section(line):
                continue
            yield line


def split_data(data, interval=100):

    splitted_data = []
    num_data = len(data)
    num_intvls = num_data/interval

    for nint in range(num_intvls):
        split_data.append( data[nint*interval:(nint+1)*interval] )

    # The remains of data are thrown away.
    return splitted_data

def cal_average_with_scalar(values_traj):

    # natom times ntraj
    natom = len(values_traj[0])
    ntraj = len(values_traj)

    sum_values = numpy.zeros([natom])
    for values in values_traj:
        sum_values += values

    return sum_values / ntraj

def cal_average_with_tensor(tensors):

    sum_tensor = numpy.zeros(tensors[0].shape)

    for t in tensors:
        sum_tensor += t

    nten = len(tensors)
    return sum_tensor/nten

def cal_rmsf(tensors, average_tensor):

    rmsf_tensor = numpy.zeros(average_tensor.shape)

    for t in tensors:
        diff_ten = average_tensor - t
        rmsf_tensor += diff_ten * diff_ten

    nten = len(tensors)

    return numpy.sqrt(rmsf_tensor) / nten

def gen_eigens(tensors):
    """Generate the calculated eigen values for each tonsor."""
    # diagonal
    for tensor in tensors:
        eigen, vec = numpy.linalg.eig(tensor)
        try:
            val = math.sqrt(eigen[0]**2 + eigen[1]**2 + eigen[2]**2)
            yield val
        except ValueError:
            print(eigen)

def get_average_after_eigen(tensors_traj):
    evalues = []
    for itraj, tensors in enumerate(tensors_traj):
        evalues.append( numpy.array(list(gen_eigens(tensors))) )
    return cal_average_with_scalar(evalues)

def get_eigen_after_average(tensors_traj):
    average_tensors = cal_average_with_tensor(tensors_traj)
    return list(gen_eigens(average_tensors))

def simplify_tensor(filename, fns, labels='', snapshot=False, **kwds):
    """Show and make figure for stress ratio."""

    if len(fns) == 0:
        pass

    elif len(labels.split(',')) != len(fns):
        print(labels.split(','))
        print(fns)
        print("The number of labels and filenames must be same.")
        quit()

    # for total stress tensor
    tot_parser = TensorParser(filename)
    names, tot_data = tot_parser.parse()
    if snapshot:
        avetot_evalues = get_average_after_eigen(tot_data)
    else:
        avetot_evalues = get_eigen_after_average(tot_data)

    # for the stress tensor for each component
    evalues_list = []
    labels = [label.strip() for label in labels.split(',') ]
    for label, fn in zip(labels, fns):
        parser = TensorParser(fn)
        names, data = parser.parse()
        if snapshot:
            evalues_list.append( get_average_after_eigen(data) )
        else:
            evalues_list.append( get_eigen_after_average(data) )

    # output the eigen value of stress tensor
    if len(evalues_list) == 0:
        header = '{id:9s}  {name:>4s} {total:>12s}'.format(
                id='#label id', name='name', total='total')
        print(header)
        fmt = '{id:>9d}  {name:<4s} {total:12.7f}'

        for name, tot_evalue in zip(names, avetot_evalues):

            # name
            atom_id, aname = name.split('_')

            # output one line data
            print(fmt.format(id=int(atom_id), name=aname, total=tot_evalue))

    else:
        header = '{id:9s}  {name:>4s} {total:>12s}'.format(
                id='#label id', name='name', total='total') + ' '.join(
                        '{:>10s}'.format(label) for label in labels)
        print(header)
        fmt = '{id:>9}  {name:<4s} {total:12.7f}' + ' '.join(
                '{'+label+':10.2f}' for label in labels)

        for tuples in zip(names, avetot_evalues, *evalues_list):

            name, tot, others = tuples[0], tuples[1], tuples[2:]

            # name
            atom_id, aname = name.split('_')

            # each evalue
            if not isinstance(others, tuple):
                others = [others]
            label_to_evalues = dict(zip(labels, others))

            # output one line data
            print(fmt.format(id=atom_id, name=aname,
                total=tot, **label_to_evalues))


if __name__ == '__main__':
    from curp.console import arg_simplify, exec_command

    parser = arg_simplify()
    exec_command(parser)

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

def cal_average(tensors):

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

def cal_contribution_ratio(vals):
    """Calculate the contribution ratio each component."""
    total = sum(vals)
    return [v/total * 100.0 for v in vals]

def gen_eigens(tensors):
    """Generate the calculated eigen values for each tonsor."""
    # diagonal
    for tensor in tensors:
        eigen, vec = numpy.linalg.eig(tensor)
        yield math.sqrt(eigen[0]**2 + eigen[1]**2 + eigen[2]**2)


if __name__ == '__main__':

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(
            'Show and make figure for stress ratio.')
        print('simplify_tensor total_file "label1,label2,label3,label4,..."',
                'label1_data,label2_data,...')

    # add argument definitions
    parser.add_argument(
            '-i', '--input-data', dest='data_filename', required=True,
            help='specify input filename for the stress data.')
    parser.add_argument(
            '-l', '--label-line', dest='label_line',
            required=False, default='',
            help='specify label string.')
    parser.add_argument(
            '-o', '--outputs', dest='fns',
            required=False, action='store_true', default=False,
            help='specify filenames.')

    # make arguments
    args = parser.parse_args()

    if len(parser.label_line.split(',')) != len(parser.fns.split(',')):
        print("The number of labels and filenames must be same.")
        quit()

    filename = args.data_filename
    label_line = args.label_line
    fns = args.split(',')

    # for total stress tonsor
    tot_parser = TensorParser(filename)
    names, tot_data = tot_parser.parse()
    tot_data_spl = splitted_data(tot_data)
    for tot_data in tot_data_spl:
        ave_tensors = cal_average(tot_data)
        tot_evalues = list(gen_eigens(ave_tensors))

    # for each component of stress tensor
    labels = [ label.strip() for label in label_line.split(',') ]
    
    evalues = []
    for lebel, fn in zip(labels, fns):
        parser = TensorParser(fn)
        names, data = parser.parse()
        average_tensors = cal_average(data)
        evalues.append( list(gen_eigens(average_tensors)) )

    # calculate for each interval
    for tot_data in -w
        ave_tensors = cal_average(tot_data)
        tot_evalues = list(gen_eigens(ave_tensors)


    # output the eigen value of stress tensor
    if len(evalues) == 0:
        header = '{id:9s}  {name:>4s} {total:>12s}'.format(
                id='#label id', name='name', total='total')
        print(header)
        fmt = '{id:>9d}  {name:<4s} {total:12.7f}'

        for name, tot_evalue in zip(names, tot_evalues):

            # name
            atom_id, aname = name.split('_')

            # output one line data
            print(fmt.format(id=int(atom_id), name=aname, total=tot_evalue))

    else:
        header = '{id:9s}  {name:>4s} {total:>12s}'.format(
                id='#label id', name='name', total='total') + ' '.join(
                        '{:>10s}'.format(label) for label in labels)
        print(header)
        fmt = '{id:>9d}  {name:<4s} {total:12.7f}' + ' '.join(
                '{'+label+':10.2f}' for label in labels)

        for tuples in zip(names, tot_evalues, *evalues):

            name, tot_evalue, evalue_tup = tuples[0], tuples[1], tuples[2:]

            # name
            atom_id, aname = name.split('_')

            # each evalue
            if not isinstance(evalue_tup, tuple):
                evalue_tup = [evalue_tup]
            ratios = cal_contribution_ratio(evalue_tup)
            label_to_ratios = dict(zip(labels, ratios))

            # output one line data
            print(fmt.format(id=int(atom_id), name=aname,
                total=tot_evalue, **label_to_ratios))


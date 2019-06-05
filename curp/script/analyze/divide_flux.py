
from __future__ import print_function
import os
import math
import numpy


import gzip
class FluxParser:

    section_end = '%end'

    def __init__(self, filename):

        if filename.endswith('.gz'):
            file = gzip.open(filename, 'rb')
        else:
            file = open(filename, 'rb')

        self.__file = file

        # make the generator
        self.__gen_snapshot_lines = self.parse_content(file)

    def __next__(self):
        try:
            lines = self.__gen_snapshot_lines.next()
            don_acc_pairs, fluxes = self.parse_onesnap(lines)
        except StopIteration:
            self.__file.close()
            raise

        return don_acc_pairs, numpy.array(fluxes)

    next = __next__ # for python 2.x

    def __iter__(self):
        return self

    def parse_header(self):
        lines = []
        for line in self.__file:
            if line[0] == '%':
                lines.append(line)
            else:
                break

        self.__file.seek(0)

        for line in lines:
            if line.startswith('%label'):
                label_line = line.replace('%label', '')

        return label_line.split()[2:]

    def close(self):
        self.__file.close()

    def parse_content(self, file):
        gen_lines = self.gen_optimized_lines(file)
        for line in gen_lines:
            if not self.is_section(line): continue

            flagname = line[1:].strip()
            if flagname != 'data': continue

            yield list( self.gen_section_lines(gen_lines) )

    def parse_onesnap(self, lines):
        """Parse the lines for one snapshot."""
        don_acc_pairs = []
        fluxes = []
        for line in lines:
            cols = line.split()
            don_acc_pairs.append( (cols[0], cols[1]) )

            values = [float(c) for c in cols[2:]]
            fluxes.append(values)

        return don_acc_pairs, numpy.array(fluxes)

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


class FluxWriter:

    def __init__(self, don_name, acc_name, base_fn, dt, labels):
        dirpath, fname = os.path.split(base_fn)
        filename = '{don}_{acc}.{base_fn}'.format(
                don=don_name, acc=acc_name, base_fn=fname)

        self.__filepath = os.path.join(dirpath, filename)
        self.__dt = dt
        self.__fmt = ('{:>12.4f}' + len(labels)*' {:16.8e}')
        self.write_header(don_name, acc_name, labels)

    def write_header(self, don_name, acc_name, labels, natom=None):
        file = self.open(self.__filepath)
        file.write('#REMARK  Flux Data for Acceptor <== Donor\n')
        file.write('#REMARK  Acceptor = {}'.format(acc_name) + '\n')
        file.write('#REMARK  Donor    = {}'.format(don_name) + '\n')
        # file.write('#REMARK  The number of atoms = {}'.format(natom) +
        info_line = '#REMARK time'
        for l in labels:
            info_line += ' {label:>16}'.format(label=l)

        file.write(info_line + '\n')
        file.close()

    def write(self, istep, flux):
        time = istep * self.__dt
        file = self.open(self.__filepath)
        file.write( self.__fmt.format(time, *flux) + '\n')
        file.close()

    def open(self, filepath):
        if filepath.endswith('.gz'):
            file = gzip.open(filepath, 'ab')
        else:
            file = open(filepath, 'ab')
        return file

def gen_flux_from_parsers(filenames):
    for fn in filenames:
        parser = FluxParser(fn)
        print('reading the data from {}'.format(fn))
        for pairs, fluxes in parser:
            yield pairs, fluxes

def gen_target_pairs(donor_line, acceptor_line, don_acc_pairs):

    if donor_line != '':
        target_dons = [ don for don in donor_line if don != '' ]
    else:
        donors = [ don for don, acc in don_acc_pairs ]
        target_dons = sorted( set(donors), key=donors.index )

    if acceptor_line != '':
        target_accs = [ acc for acc in acceptor_line if acc != '' ]
    else:
        acceptors = [ acc for don, acc in don_acc_pairs ]
        target_accs = sorted( set(acceptors), key=acceptors.index )

    print('target_dons', target_dons)
    print('target_accs', target_accs)

    for don, acc in don_acc_pairs:
        if don not in target_dons and acc not in target_accs: continue
        yield don, acc

def write_summary(filename, don_acc_pairs, labels, fluxes, wtype):

    basename, ext = os.path.splitext(filename)
    filepath = basename + '_' + wtype + ext

    open_ = gzip.open if filepath.endswith('.gz') else open

    # write header
    with open_(filepath, 'wb') as file:
        print('writing the {} data into {}'.format(wtype, filepath))
        file.write('%title {} flux\n'.format(wtype))
        label_line = '% {:>10} {:>12}'.format('donor', 'acceptor')
        for l in labels:
            label_line += ' {:>12s}'.format(l)
        file.write(label_line + '\n')

        fmt = '{:>12} {:>12}' + ' {:12.7f}'*len(labels) + '\n'
        for (don, acc), flux in zip(don_acc_pairs, fluxes):
            file.write( fmt.format(don, acc, *flux))

def write_avg_rms(filename, don_acc_pairs, avg_fluxes, rms_fluxes):

    basename, ext = os.path.splitext(filename)
    filepath = basename + '_avg+rms' + ext

    open_ = gzip.open if filepath.endswith('.gz') else open

    # write header
    with open_(filepath, 'wb') as file:
        print('writing the {} data into {}'.format('avg+rms', filepath))
        file.write('%title average and rms of the flux\n')
        label_line = '% {:>10} {:>12} {:>12} {:>12}'.format(
                'donor', 'acceptor', 'average', 'rms')
        file.write(label_line + '\n')

        fmt = '{:>12} {:>12} {:12.7f} {:12.7f}' + '\n'
        for (don, acc), avg, rms in zip(don_acc_pairs, avg_fluxes, rms_fluxes):
            file.write( fmt.format(don, acc, avg[0], rms[0]))


def divide_flux(flux_fns, output_fn, dt=1.0, donor_line='', acceptor_line='',
         column_line='', **kwds):
    """Divide flux file.

    Divide the flux data of all into the flux data for each donor and
    acceptor.
    """

    # prepare
    one_parser = FluxParser(flux_fns[0])
    all_labels = one_parser.parse_header()
    don_acc_pairs, fluxes = one_parser.next()
    one_parser.close()
    nvalue = len(all_labels)

    # parse colums list. For example, '1,4,6,9'
    if column_line == '':
        icolums = range(1, nvalue+1)
    else:
        icolums = [ int(col) for col in column_line.split(',') ]

    print('icolums', icolums)

    labels = []
    for icol in icolums:
        label = all_labels[icol-1]
        labels.append( label )

    print('labels', labels)

    # create average fluxes and fluxes**2 array
    npair = len(don_acc_pairs)
    print('npair', npair)
    fluxes_sum  = numpy.zeros([npair, nvalue], numpy.float)
    fluxes2_sum = numpy.zeros([npair, nvalue], numpy.float)

    # get target pairs
    target_pairs = list(gen_target_pairs(donor_line, acceptor_line,
                                          don_acc_pairs))
    target_indexes = []
    for index, pair in enumerate(target_pairs):
        target_indexes.append(index)

    print('target_pairs', target_pairs)
    print('target_indexes', target_indexes)

    # make writers
    ipair_to_writer = []
    for don_name, acc_name in target_pairs:
        writer = FluxWriter(don_name, acc_name,
                            output_fn, float(dt), labels)
        ipair_to_writer.append( writer )

    # parse and write the energy flux
    parsers = gen_flux_from_parsers(flux_fns)
    for istep, (pairs, fluxes) in enumerate(parsers):
        fluxes_sum  += fluxes
        fluxes2_sum += fluxes*fluxes

        for index, writer in zip( target_indexes, ipair_to_writer ):
            writer.write( istep, fluxes[index] )

    # postprocessing
    # for writer in ipair_to_writer:
    #     writer.close()

    avg_fluxes = fluxes_sum / float(istep+1)
    avg_fluxes2 = fluxes2_sum / float(istep+1)
    rms_fluxes = numpy.sqrt( avg_fluxes2 - avg_fluxes**2 )

    if len(labels) == 1:
        write_avg_rms(output_fn, target_pairs, avg_fluxes, rms_fluxes)
    else:
        # calculate average fluxes
        write_summary(output_fn, target_pairs, labels, avg_fluxes, 'avg')
        write_summary(output_fn, target_pairs, labels, rms_fluxes, 'rms')

    print()
    print(80*'-')
    print('    Script that divide flux data finished completely.')
    print(80*'-')


if __name__ == '__main__':
    from curp.console import arg_divide_flux, exec_command

    parser = arg_divide_flux()
    exec_command(parser)

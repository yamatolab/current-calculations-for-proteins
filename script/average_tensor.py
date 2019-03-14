
import sys
import numpy

section_end = '%end'

def parse_content(file):

    tensors = []

    gen_lines = gen_optimized_lines(file)
    sections = []
    for line in gen_lines:
        if not is_section(line): continue

        flagname = line[1:].strip()
        if flagname == 'data':
            lines = list( gen_section_lines(gen_lines) )
            names, tensor = parse_data(lines)
            tensors.append( numpy.array(tensor) )

    return names, tensors

def parse_data(lines):
    names = []
    tensor = []
    for line in lines:
        cols = line.split()
        names.append(cols[0])
        tensor.append( [float(c) for c in cols[1:]] )

    return names, tensor

def gen_optimized_lines(file):
    # first flag
    for line in file:
        if is_section(line):
            yield line
            break
        yield line

    for line in file:
        if is_section(line):
            yield section_end
            yield line
        yield line

    else:
        yield section_end

def split_content(gen):
    lines = []
    for line in gen:
        if line.startswith('%'):
            pass

def is_section(line):
    if line == section_end:
        return False
    elif line.startswith('%'):
        return True
    else:
        return False

def gen_section_lines(gen):
    for line in gen:
        if line == section_end:
            break
        elif is_section(line):
            continue
        yield line



            


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

if __name__ == '__main__':

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser('Get the averaged tensor.')

    # add argument definitions
    parser.add_argument(
            '-i', '--input-data', dest='data_filename', required=True,
            help='specify input filename for the stress data.')

    # make arguments
    args = parser.parse_args()

    filename = args.data_filename
    with open(filename, 'rb') as file:
        names, tensors = parse_content(file)

    ave_tensor = cal_average(tensors)
    rmsf_tensor = cal_rmsf(tensors, ave_tensor)
    for name, avet, rmsft in zip(names, ave_tensor, rmsf_tensor):
        print(name)
        print(avet)
        print(rmsft)

    



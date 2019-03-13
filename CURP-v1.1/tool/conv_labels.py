#! /usr/bin/env python2
from __future__ import print_function
import sys, os

class LabelConverterParser:

    def __init__(self, filename):
        self.parse(filename)

    def parse(self, filename):
        file = open(filename, 'rb')
        label_to_newlabel = {}
        for line in file:
            line = line.strip()
            if not line: continue
            if line.startswith('#'): continue

            key, newlabel = line.split(':')
            label_to_newlabel[key.strip()] = newlabel.strip()

        self.__label_to_newlabel = label_to_newlabel
        file.close()

    def get_label_to_newlabel(self):
        return self.__label_to_newlabel


if __name__ == '__main__':

    ec_file  = sys.stdin
    label_fp = sys.argv[1]

    if not os.path.exists(label_fp):
        for line in ec_file:
            print(line.strip())
        exit()

    label_parser = LabelConverterParser(label_fp)
    label_to_newlabel = label_parser.get_label_to_newlabel()

    for line in ec_file:
        line = line.strip()
        if not line: continue

        if line.startswith('#'):
            print(line.strip())
            continue

        cols = line.split()
        don_line, acc_line, rest = cols[0], cols[1], cols[2:]

        don_line_ = label_to_newlabel.get(don_line)
        if don_line_: don_line = don_line_

        acc_line_ = label_to_newlabel.get(acc_line)
        if acc_line_: acc_line = acc_line_

        print(don_line, acc_line, ' '.join(rest))

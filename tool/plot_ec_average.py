#! /usr/bin/env python2
from __future__ import print_function
import os, sys
import numpy
import pylab as pl

def load_data(filename):
    with open(filename, 'rb') as file:

        don_acc_lines = []
        ecs = []

        lines = (line for line in file
                if not line.startswith('#')
                if not line.isspace() )

        for line in lines:
            cols = line.split()
            don_acc_lines.append( (cols[0], cols[1]) )
            ecs.append( float(cols[2]) )

    return don_acc_lines, ecs

if __name__ == '__main__':

    base_fn = sys.argv[1]
    rel_fns = sys.argv[2:]

    pl.figure(figsize=(10,6))


    base_lines, base_ecs = load_data(base_fn)
    base_ecs = numpy.array(base_ecs)

    markers = ['^', 'D', 'o']
    colors  = ['red', 'green', 'blue']

    for fn, color, marker in zip(rel_fns, colors, markers):
        rel_lines, rel_ecs = load_data(fn)
        rel_ecs = numpy.array(rel_ecs)

        xs = rel_ecs/base_ecs

        pl.plot(xs, alpha=0.7, linestyle='_',
                color=color, marker=marker, markersize=4)

        for i, (x, base_ec, rel_ec) in enumerate(zip(xs, base_ecs, rel_ecs)):
            if x >= 1.2 or x <= 0.7:
                print(base_lines[i], x, base_ec, rel_ec)


    pl.ylim((0,3.0))
    pl.ylabel(r'$\langle irEC_{ij} \rangle_x / \langle irEC_{ij} \rangle_{100} $')
    pl.xlabel(r'Residue pair number')
    # pl.legend(('x=10', 'x=20', 'x=50'))

    rel_fn = rel_fns[-1]
    name = os.path.splitext(os.path.basename(rel_fn))[0]

    pl.savefig('avg-'+name+'.png', format='png', dpi=300)
    
    # pl.show()

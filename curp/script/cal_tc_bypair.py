#! /usr/bin/env python
from __future__ import print_function
import os, sys
import math
import numpy
import gzip

# fortran module
from curp.script.lib_acf import cal_acf

# Exception
class NumFrameError(Exception): pass

def parse_options():
    """Parse and get the command line options."""

    # make argument parser
    import argparse
    parser = argparse.ArgumentParser(
        description= ('Calculate auto-correlation function data '
            'from energy flux datas.') )

    # add argument definitions
    parser.add_argument(
            'flux_fns', nargs='+',
            help=('specify filenames of flux data. '
                'ex.) flux.dat0001, flux.dat0002, ...') )

    parser.add_argument(
            '-a', '--output-acf', dest='acf_fn', required=False, default='',
            help='specify a filename to output acf.')

    parser.add_argument(
            '-D', '--without-dispersion', dest='enable_disp',
            action='store_false', required=False,
            help='calculate dispersion for tc.')

    parser.add_argument(
            '-t', '--output-tcs', dest='tcs_fn', required=False, default='',
            help=('specify a filename to output'
                  'time series of transport coefficient.'))

    # parser.add_argument(
            # '-p', '--with-precise', dest='enable_precise',
            # action='store_true', required=False,
            # help='calculate dispersion in precise.')

    parser.add_argument(
            '-f', '--frame-range', dest='frame_range',
            type=int, nargs=3, required=True,
            help='frame range; first_frame, last_frame, intervalstep.')

    parser.add_argument(
            '-s', '--shift', dest='frm_shift',
            type=int, required=False, default=1,
            help='shift frame to collect samples.')

    parser.add_argument(
            '--sample-number', dest='nsample',
            type=int, required=False, default=0,
            help=('The number of samples for one flux data file. '
                'Default value is 0 that present to make samples '
                'as much as possible'))

    parser.add_argument(
            '-n', '--with-norm', dest='enable_norm',
            action='store_true', required=False,
            help='flag to normalize acf.')

    parser.add_argument(
            '-dt', '--dt', dest='dt',
            type=float, required=False, default=1.0,
            help='t of between the neighbour frames. The unit is in ps')

    parser.add_argument(
            '-c', '--coefficient', dest='coef',
            type=float, required=False, default=1.0,
            help='write coefficient times acf.')

    parser.add_argument(
            '-o', '--output-log', dest='log_fn', required=False,
            default='tc.log',
            help='specify a filename for log.')

    parser.add_argument(
            '-v', '--vervose', dest='use_debug',
            action='store_true', required=False,
            help='turn on debug mode.')

    # make arguments
    options = parser.parse_args()

    return options

def load_fluxdata(flux_fn):

    if flux_fn.endswith('.gz'):
        file = gzip.open(flux_fn, 'rb')
    else:
        file = open(flux_fn, 'rb')

    ifrm = 0
    xss = []
    for line in file:
        if line.startswith('#'): continue

        cols = line.split()
        ifrm += 1
        xss.append( float(cols[1]) )

    file.close()
    return numpy.array(xss)

def cal_tc(acf, dt, coef=1.0):
    """
    Calculate final transport  coefficient.
    old version: return coef * dt*1000.0 * numpy.sum(acf)
    """
    return coef * dt*1000.0 * numpy.trapz(acf,axis=1)[0]

def cal_tcs(acf, dt, coef=1.0):
    """Calculate time series of transport coefficient."""
    return coef * dt*1000.0 * numpy.add.accumulate(acf)

def write_acf(times, acf, opts):
    """Write auto-correlation function data to given file."""

    # open file
    acf_fn = opts.acf_fn
    if acf_fn.endswith('.gz'):
        file = gzip.open(acf_fn, 'wb')
    else:
        file = open(acf_fn, 'wb')

    # get donor and acceptor id and name
    flux_fn = opts.flux_fns[0]
    don_acc_line = os.path.basename(flux_fn).split('.')[0]
    don_id, don_name, acc_id, acc_name = don_acc_line.split('_')
    don_id, acc_id = int(don_id), int(acc_id)

    # write header
    file.write('#donor    {}{}\n'.format(don_name, don_id))
    file.write('#acceptor {}{}\n'.format(acc_name, acc_id))
    file.write('#label {:>8}  {:>22}\n'.format('time(ps)','acf(kcal/mol/fs)^2'))

    # write data
    for time, a in zip(times, acf[:, 0]):
        file.write("{:15.4f}  {:22.14e}\n".format(time, a))

    # close acf data file
    file.close()

def write_tc(tc, tc_disp, don_acc_line):
    """Write transport coefficient and dispersion to standard output."""
    if tc_disp:
        print( don_acc_line, tc, '+-', tc_disp)
    else:
        print( don_acc_line, tc)

def get_don_acc_line(flux_fn):
    don_acc_line = os.path.basename(flux_fn).split('.')[0]
    don_id, don_name, acc_id, acc_name = don_acc_line.split('_')
    return "{}_{} {}_{}".format( don_id, don_name, acc_id, acc_name)

def write_tc_traj(itrj, tc, don_acc_line):
    """Write transport coefficient and dispersion to standard output."""


def write_acf(times, acf, acf_fn, flux_fn):

    # open file
    if acf_fn.endswith('.gz'):
        file = gzip.open(acf_fn, 'wb')
    else:
        file = open(acf_fn, 'wb')

    # get donor and acceptor id and name
    don_acc_line = os.path.basename(flux_fn).split('.')[0]
    don_id, don_name, acc_id, acc_name = don_acc_line.split('_')
    don_id, acc_id = int(don_id), int(acc_id)

    # write header
    file.write('#donor    {}{}\n'.format(don_name, don_id))
    file.write('#acceptor {}{}\n'.format(acc_name, acc_id))
    file.write('#label {:>8}  {:>22}\n'.format(
        'time(ps)', 'acf[(kcal/mol/fs)^2]'))

    # write data
    for time, a in zip(times, acf[:, 0]):
        file.write("{:15.4f}  {:22.14e}\n".format(time, a))

    # close tcs data file
    file.close()

def write_tcs(times, tcs, disps, tcs_fn, flux_fn):

    # open file
    if tcs_fn.endswith('.gz'):
        file = gzip.open(tcs_fn, 'wb')
    else:
        file = open(tcs_fn, 'wb')

    # get donor and acceptor id and name
    don_acc_line = os.path.basename(flux_fn).split('.')[0]
    don_id, don_name, acc_id, acc_name = don_acc_line.split('_')
    don_id, acc_id = int(don_id), int(acc_id)

    # write header
    file.write('#donor    {}{}\n'.format(don_name, don_id))
    file.write('#acceptor {}{}\n'.format(acc_name, acc_id))
    file.write('#label {:>8}  {:>22} {:>22}\n'.format(
        'time(ps)', 'tcs[(kcal/mol)^2/fs])', 'dispersion' ) )

    # write data
    for time, tc, disp in zip(times, tcs[:,0], disps[:,0]):
        file.write("{:15.4f}  {:22.14e} {:22.14e}\n".format(time, tc, disp))

    # close tcs data file
    file.close()

class Log:
    def __init__(self, filename):
        self.__fn = filename

        # initialize log file
        with open(filename, 'wb') as file:
            pass

    def write(self, line):
        with open(self.__fn, 'ab') as file:
            file.write(str(line) + '\n')

def sum_gen(gen):
    first = gen.next()
    for g in gen:
        first += g

    return first


def main():

    # parse command line options
    opts = parse_options()

    log = Log(opts.log_fn)
    if opts.use_debug:
        log.write(opts)

    frm_first, frm_last, frm_interval = opts.frame_range
    ncom = 1

    # prepare
    nfile = len(opts.flux_fns)
    nacf = (frm_last - frm_first)/frm_interval + 1

    # generator acf by data file
    def gen_acf(flux_fns, tag=''):
        for ifn_1, flux_fn in enumerate(flux_fns):
            fluxes = load_fluxdata(flux_fn)

            if len(fluxes) < nacf:
                raise NumFrameError(("\n"+
                    "   Number of flux data is {},\n" +
                    "   but given frame range is {}.")
                    .format(len(fluxes), nacf))

            acf = cal_acf(fluxes,nacf,frm_first,frm_last,frm_interval,
                    opts.frm_shift, norm=opts.enable_norm, nsample=opts.nsample)
            log.write('[{tag}] finished {number}'.
                    format(tag=tag, number=ifn_1+1))
            yield acf

    # set coefficient
    dt = opts.dt * opts.frame_range[2] # in ps
    times = numpy.arange(0.0, nacf*dt, dt)

    ################################
    # final tc and dispersion only #
    ################################
    don_acc_line = get_don_acc_line(opts.flux_fns[0])
    if (not opts.acf_fn) and opts.enable_disp:

        acf_iter = gen_acf(opts.flux_fns, 'tc+disp')
        tc_ary = numpy.array( [cal_tc(acf, dt, opts.coef) for acf in acf_iter] )
        for itraj, tc in enumerate(tc_ary): # by trajectory(file)
            print("%itraj="+str(itraj), don_acc_line, tc)

        # calculate tc and dispersion
        tc_avg  = numpy.sum(tc_ary) / nfile
        tc2_avg = numpy.sum(tc_ary*tc_ary) / nfile

        # average tc from tc array
        tc = tc_avg
        tc_disp = math.sqrt( tc2_avg - tc_avg**2 )
        print( don_acc_line, tc, '  ', tc_disp)

    ####################
    # acf and final tc #
    ####################
    elif opts.acf_fn and (not opts.enable_disp):
        acf_iter = gen_acf(opts.flux_fns, 'acf+tc')

        # get average acf and tc
        acf_avg = sum(acf_iter) / nfile
        # calculate tc from averaged acf
        tc = cal_tc(acf_avg, dt, opts.coef)

        write_acf(times, acf_avg, opts.acf_fn, opts.flux_fns[0])
        print( don_acc_line, tc)

    ##################################
    # acf and final tc +- dispersion #
    ##################################
    elif opts.acf_fn and opts.enable_disp:
        acf_iter = gen_acf(opts.flux_fns, 'acf+tc+disp')

        tc_list = []
        acf_sum = numpy.zeros([nacf,ncom])
        for acf in acf_iter:
            acf_sum += acf
            tc_list.append(cal_tc(acf, dt, opts.coef))

        # calculate tc and dispersion
        tc_ary  = numpy.array(tc_list)
        for itraj, tc in enumerate(tc_ary):
            print("%itraj="+str(itraj), don_acc_line, tc)

        tc_avg  = numpy.sum(tc_ary) / nfile
        tc2_avg = numpy.sum(tc_ary*tc_ary) / nfile

        # get average_acf and tc dispersion
        acf_avg = acf_sum / nfile
        tc      = cal_tc(acf_avg, dt, opts.coef)
        # tc_avg and tc muse be equal !
        tc_disp = math.sqrt( tc2_avg - tc_avg**2 )

        write_acf(times, acf, opts.acf_fn, opts.flux_fns[0])
        print( don_acc_line, tc, '  ', tc_disp)

    #################################################
    #  get tcs( time-series of tc ) and dispersions #
    #################################################
    if opts.tcs_fn: # and (not opts.enable_precise):
        acf_iter = gen_acf(opts.flux_fns, 'tcs')

        tcs_sum  = numpy.zeros([nacf,ncom])
        tcs_sum2 = numpy.zeros([nacf,ncom])
        for acf in acf_iter:
            tcs = cal_tcs(acf, dt, opts.coef)
            tcs_sum  += tcs
            tcs_sum2 += tcs * tcs

        tcs_avg  = tcs_sum / nfile
        tcs2_avg = tcs_sum2 / nfile

        # dispersion of tc based on <x^2> - <x>^2
        tcs_disp = numpy.sqrt( tcs2_avg - tcs_avg**2 )
        # tc_disp and tcs_disp[-1] must be equal !
        if opts.use_debug: log.write('tc = {:22.14e} {:22.14e}'
            .format(tc, tcs_avg[-1][0]))
        if opts.use_debug: log.write('disp = {:22.14e} {:22.14e}'
            .format(tc_disp, tcs_disp[-1][0]))

        write_tcs(times, tcs_avg, tcs_disp, opts.tcs_fn, opts.flux_fns[0])

    ############################################################
    #  get tcs( time-series of tc ) and dispersions in precise #
    ############################################################
    # if opts.tcs_fn and opts.enable_precise:
        # tcs_avg = cal_tcs(acf_avg, dt, opts.coef)

        # acf_iter = gen_acf(opts.flux_fns)

        # tcs2_sum  = numpy.zeros([nacf,ncom])

        # def gen_diff2(tcs_avg, acf_iter):
            # for acf in acf_iter:
                # tcs = cal_tcs(acf, dt)
                # yield (tcs - tcs_avg)**2

        # # dispersion of tc based on <(dx)^2>
        # tcs_disp = numpy.sqrt( sum(gen_diff2(tcs_avg, acf_iter)) / nfile )

        # write_tcs(times, tcs_avg, tcs_disp, opts.tcs_fn, opts.flux_fns[0])

    # write
    # write_tc(acf, opts)

if __name__ == '__main__':
    main()

"""clog - Defines functions of what to print in the logs.

If log_level is set to 1 (WARN):
  only warning messages will be printed.
If log level is set to 2 (default) (INFO):
  warning and info messages will be printed.
If log level is set to 4 (DEBUG):
  warning, info and debug messages will be printed.
"""

from __future__ import print_function


WARN = 1
INFO = 2
DEBUG = 4 

print_function = print
log_level = INFO

title_fmt = '\n' + 80*'-' + '\n' + 8*' '+ '{}' '\n' + 80*'-' + '\n'
title_fmt_info = '\n' + 80*'-' + '\n' + '{:>6}. {}' '\n' + 80*'-' + '\n'
count = 0
log_frequency = 1
cur_istep = 0

def set_log_frequency(frequency):
    global log_frequency
    log_frequency = frequency

def set_curstep(istp):
    global cur_istep 
    cur_istep = istp

def print(*args, **kwds):
    print_function(*args, **kwds)

#Debug functions

def debug(*args, **kwds):
    if log_level >= DEBUG: # log_level = 4
        print_function(*args, **kwds)

def debug_cycle(*args, **kwds):
    """Log informations depending on log_frequency."""
    if cur_istep in [0,1] or cur_istep%log_frequency==0:
        debug(*args, **kwds)

def debug_title(title):
    debug(title_fmt.format(title))

def is_debug():
    return log_level == DEBUG

#Info functions

def info(*args, **kwds):
    if log_level >= INFO:
        print_function(*args, **kwds)

def info_cycle(*args, **kwds):
    """Log informations depending on log_frequency."""
    if cur_istep in [0,1] or cur_istep%log_frequency==0:
        info(*args, **kwds)

def info_title(title):
    global count
    count += 1
    info(title_fmt_info.format(count, title))

def is_info():
    return log_level == INFO

#Warn functions

def warn(*args, **kwds):
    if log_level >= WARN:
        print_function(*args, **kwds)

def warn_title(title):
    warn(title_fmt.format(title))
def is_warn():
    return log_level == WARN

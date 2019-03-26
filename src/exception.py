#  -*- encoding: utf-8 -*-
"""exception - Exception handler"""
import os, sys
import traceback

from datetime import datetime
from pprint import pformat

raise_exception = True

class CurpException(Exception):

    def __init__(self, data=None, log_fn='log.txt', tb=False):
        BaseException.__init__(self)
        self.message = 'Curp Exception'
        self.log_fn = log_fn
        self.date = datetime.now()
        self.__data = data
        self.__traceback = tb

#    def __repr__(self):
#        return 'CurpException'

    def __str__(self):
        mes = pformat(self.__data) if self.__data else ''
        if self.__traceback:
            mes += '\n'
            mes += traceback.format_exc()

        return mes

    def log(self):
        print self.__data

    def logerror(self):
        logfile = open(self.log_fn, 'a')
        logfile.write(str(self.date) + ' : ' + str(self))
        logfile.close()

    def getMessage(self):
        return self.__data


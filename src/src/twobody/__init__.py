import os, sys

topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)

from exception import CurpException
class PotentialNotFoundError(CurpException): pass

################################################################################
def get_calculator(forcefield):
    if forcefield in ('amber-base','amber94','amber96','amber99',
            'amber99SB','amber03', 'amber10', 'amber12SB'):
        import amberbase
        return amberbase.TwoBodyForce

    else:
        raise PotentialNotFoundError(forcefield)


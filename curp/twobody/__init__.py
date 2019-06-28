from curp.exception import CurpException


class PotentialNotFoundError(CurpException): pass


def get_calculator(forcefield):
    if forcefield in ('amber-base', 'amber94', 'amber96', 'amber99',
                      'amber99SB', 'amber03', 'amber10', 'amber12SB'):
        from curp.forcefield import amberbase
        return amberbase.TwoBodyForce

    else:
        raise PotentialNotFoundError(forcefield)

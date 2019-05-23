import current
import flux

def get_calculator(setting):
    """Return the calculator object by setting."""

    # get volume setting
    method = setting.curp.method

    if method == 'momentum-current':
        obj = current.StressCurrentCalculator
    elif method == 'stress-flux':
        obj = flux.StressFluxCalculator
    elif method == 'energy-current':
        obj = current.EnergyCurrentCalculator
    elif method == 'energy-flux':
        obj = flux.EnergyFluxCalculator
    elif method == 'heat-flux':
        obj = flux.HeatFluxCalculator
    else:
        raise Exception

    return obj



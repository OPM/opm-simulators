import os
from contextlib import contextmanager

@contextmanager
def pushd(path):
    cwd = os.getcwd()
    if not os.path.isdir(path):
        os.makedirs(path)
    os.chdir(path)
    yield
    os.chdir(cwd)

ENABLE_ASYNC_ECL_OUTPUT_FLAG = '--enable-async-ecl-output=false'

def create_black_oil_simulator(*args, **kwargs):
    """Create BlackOilSimulator with test-safe default arguments.

    Automatically disables async ECL output to prevent race conditions
    with pushd context manager in tests.
    """
    from opm.simulators import BlackOilSimulator

    flag_to_add = ENABLE_ASYNC_ECL_OUTPUT_FLAG
    # Handle different constructor patterns
    if 'args' in kwargs:
        # Add our flag to existing args
        if flag_to_add not in kwargs['args']:
            kwargs['args'].append(flag_to_add)
    elif len(args) >= 4:
        # Constructor with deck, state, schedule, summary_config
        # Add args parameter
        if len(args) == 4:
            args = args + ([flag_to_add],)
        elif len(args) == 5:
            # Args already provided, add our flag
            args_list = list(args[4]) if args[4] else []
            if flag_to_add not in args_list:
                args_list.append(flag_to_add)
            args = args[:4] + (args_list,)
    else:
        # Constructor with filename - add args parameter
        kwargs['args'] = kwargs.get('args', [])
        if flag_to_add not in kwargs['args']:
            kwargs['args'].append(flag_to_add)

    return BlackOilSimulator(*args, **kwargs)

def create_gas_water_simulator(*args, **kwargs):
    """Create GasWaterSimulator with test-safe default arguments.

    Automatically disables async ECL output to prevent race conditions
    with pushd context manager in tests.
    """
    from opm.simulators import GasWaterSimulator

    flag_to_add = ENABLE_ASYNC_ECL_OUTPUT_FLAG
    # Add our flag to args
    kwargs['args'] = kwargs.get('args', [])
    if flag_to_add not in kwargs['args']:
        kwargs['args'].append(flag_to_add)

    return GasWaterSimulator(*args, **kwargs)

def create_onephase_simulator(*args, **kwargs):
    """Create OnePhaseSimulator with test-safe default arguments.

    Automatically disables async ECL output to prevent race conditions
    with pushd context manager in tests.
    """
    from opm.simulators import OnePhaseSimulator

    flag_to_add = ENABLE_ASYNC_ECL_OUTPUT_FLAG
    # Add our flag to args
    kwargs['args'] = kwargs.get('args', [])
    if flag_to_add not in kwargs['args']:
        kwargs['args'].append(flag_to_add)

    return OnePhaseSimulator(*args, **kwargs)

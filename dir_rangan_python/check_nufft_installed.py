from torch import device, cuda
from importlib.util import find_spec
from get_device import get device

def check_nufft_installed(dev: device) -> None:
    """ (from Jeff's device_handling: 20251003)
    Check for the availability of NUFFT libraries.
    """
    if dev.type == 'cuda':
        spec = find_spec('cufinufft')
        if spec is None:
            raise Exception("CUDA is requested, but cufinufft is not installed.")
    else:
        spec = find_spec('finufft')
        if spec is None:
            raise Exception("CPU is requested, but finufft is not installed.")
    pass
